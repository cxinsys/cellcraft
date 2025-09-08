"""
Snakemake Native DAG Parser

Snakemake의 네이티브 기능들(--dag, --rulegraph, --filegraph, --d3dag)을 활용하여
정확하고 효율적인 DAG 구조를 생성하는 파서 모듈.

기존의 복잡한 정규표현식 기반 파싱을 대체하여 다음과 같은 이점을 제공:
- 정확한 의존성 관계 (Snakemake 엔진 기반)
- 유지보수성 향상 (Snakemake 업데이트 자동 반영)
- 성능 최적화 (네이티브 명령어 효율성)
- 코드 복잡도 대폭 감소
"""

import os
import re
import json
import time
import logging
import subprocess
import hashlib
from typing import Dict, List, Any, Optional, Tuple, Union
from pathlib import Path
from dataclasses import dataclass
# Threading and caching removed for simplification
from datetime import datetime

# 로거 설정
logger = logging.getLogger(__name__)

@dataclass
class SnakemakeNativeError(Exception):
    """Snakemake 네이티브 파싱 관련 커스텀 예외 클래스"""
    message: str
    command: Optional[str] = None
    stderr: Optional[str] = None
    returncode: Optional[int] = None
    details: Optional[str] = None

class SnakemakeCommandRunner:
    """
    Snakemake 명령어를 안전하게 실행하는 클래스
    Docker 컨테이너 환경과 일반 환경 모두 지원
    """
    
    def __init__(self, timeout: int = 30, cwd: Optional[str] = None):
        self.timeout = timeout
        self.cwd = cwd or os.getcwd()
    
    def run_snakemake_command(self, snakefile_path: str, command_type: str, 
                            additional_args: Optional[List[str]] = None) -> Tuple[str, str, int]:
        """
        Snakemake 명령어를 실행하고 결과를 반환
        
        Args:
            snakefile_path: Snakefile 경로
            command_type: 명령어 타입 ('dag', 'rulegraph', 'filegraph', 'd3dag')
            additional_args: 추가 인수들
            
        Returns:
            (stdout, stderr, returncode) 튜플
        """
        with self._lock:
            try:
                # 기본 명령어 구성
                cmd = ['snakemake', '--snakefile', snakefile_path]
                
                # 명령어 타입에 따른 옵션 추가
                if command_type == 'dag':
                    cmd.append('--dag')
                elif command_type == 'rulegraph':
                    cmd.append('--rulegraph')
                elif command_type == 'filegraph':
                    cmd.append('--filegraph')
                elif command_type == 'd3dag':
                    cmd.append('--d3dag')
                else:
                    raise SnakemakeNativeError(f"Unsupported command type: {command_type}")
                
                # 추가 인수들 추가
                if additional_args:
                    cmd.extend(additional_args)
                
                # Dry-run으로 실행하여 안전성 확보
                cmd.extend(['--dry-run', '--quiet'])
                
                logger.info(f"Executing Snakemake command: {' '.join(cmd)}")
                
                # 명령어 실행
                result = subprocess.run(
                    cmd,
                    cwd=self.cwd,
                    capture_output=True,
                    text=True,
                    timeout=self.timeout,
                    env=self._prepare_environment()
                )
                
                logger.debug(f"Command completed with return code: {result.returncode}")
                
                return result.stdout, result.stderr, result.returncode
                
            except subprocess.TimeoutExpired:
                error_msg = f"Snakemake command timed out after {self.timeout}s"
                logger.error(error_msg)
                raise SnakemakeNativeError(
                    error_msg,
                    command=' '.join(cmd),
                    details="Command execution exceeded timeout limit"
                )
            except FileNotFoundError:
                error_msg = "Snakemake executable not found"
                logger.error(error_msg)
                raise SnakemakeNativeError(
                    error_msg,
                    command=' '.join(cmd),
                    details="Snakemake is not installed or not in PATH"
                )
            except Exception as e:
                error_msg = f"Failed to execute Snakemake command: {e}"
                logger.error(error_msg)
                raise SnakemakeNativeError(
                    error_msg,
                    command=' '.join(cmd),
                    details=str(e)
                )
    
    def _prepare_environment(self) -> Dict[str, str]:
        """실행 환경 준비"""
        env = os.environ.copy()
        
        # Python 경로 설정 (필요한 경우)
        if 'PYTHONPATH' not in env:
            env['PYTHONPATH'] = os.getcwd()
        
        # 로그 레벨 설정
        env['SNAKEMAKE_LOG_LEVEL'] = 'ERROR'
        
        return env

class DOTParser:
    """DOT 형식의 그래프 데이터를 파싱하는 클래스"""
    
    @staticmethod
    def parse_dot_graph(dot_content: str) -> Dict[str, Any]:
        """
        DOT 형식 그래프를 파싱하여 노드와 엣지 정보 추출
        
        Args:
            dot_content: DOT 형식의 그래프 내용
            
        Returns:
            노드와 엣지 정보를 포함한 딕셔너리
        """
        try:
            nodes = []
            edges = []
            
            # DOT 파일에서 노드와 엣지 추출
            lines = dot_content.strip().split('\n')
            
            for line in lines:
                line = line.strip()
                
                # 빈 줄이나 주석, DOT 키워드 무시
                if (not line or line.startswith('//') or line.startswith('digraph') or 
                    line.startswith('}') or line.startswith('graph[') or 
                    line.startswith('node[') or line.startswith('edge[')):
                    continue
                
                # 노드 정의 파싱 (예: 0[label = "prepare_data", color = "0.33 0.6 0.85", style="rounded"];)
                node_match = re.match(r'^\s*(\d+)\[([^\]]+)\]\s*;?$', line)
                if node_match:
                    node_id = node_match.group(1)
                    attributes = node_match.group(2)
                    
                    # label 속성 추출 (형식: label = "rule_name")
                    label_match = re.search(r'label\s*=\s*"([^"]+)"', attributes)
                    label = label_match.group(1) if label_match else node_id
                    
                    # 색상이나 모양 정보 추출 (있는 경우)
                    color_match = re.search(r'color\s*=\s*"([^"]+)"', attributes)
                    color = color_match.group(1) if color_match else None
                    
                    shape_match = re.search(r'shape\s*=\s*"([^"]+)"', attributes)
                    shape = shape_match.group(1) if shape_match else None
                    
                    style_match = re.search(r'style\s*=\s*"([^"]+)"', attributes)
                    style = style_match.group(1) if style_match else None
                    
                    nodes.append({
                        'id': node_id,
                        'label': DOTParser._format_label(label),
                        'shortLabel': DOTParser._create_short_label(label),
                        'type': DOTParser._infer_node_type(label),
                        'color': color,
                        'shape': shape,
                        'style': style,
                        'description': DOTParser._create_description(label)
                    })
                
                # 엣지 정의 파싱 (예: 0 -> 1 [color=grey, penwidth=2];)
                edge_match = re.match(r'^\s*(\d+)\s*->\s*(\d+)', line)
                if edge_match:
                    source = edge_match.group(1)
                    target = edge_match.group(2)
                    
                    edges.append({
                        'source': source,
                        'target': target,
                        'label': f'{source} → {target}'
                    })
            
            return {
                'nodes': nodes,
                'edges': edges,
                'format': 'DOT'
            }
            
        except Exception as e:
            logger.error(f"Error parsing DOT graph: {e}")
            raise SnakemakeNativeError(
                "Failed to parse DOT format graph",
                details=str(e)
            )
    
    @staticmethod
    def _format_label(label: str) -> str:
        """라벨을 사용자 친화적으로 포맷팅"""
        # Snakemake 규칙명에서 불필요한 부분 제거
        formatted = re.sub(r'^\w+__', '', label)  # 접두사 제거
        formatted = formatted.replace('_', ' ').title()
        return formatted
    
    @staticmethod
    def _create_short_label(label: str) -> str:
        """사각형 노드에 표시할 짧은 라벨 생성"""
        # 긴 라벨을 줄바꿈으로 분할
        words = label.replace('_', ' ').split()
        if len(words) > 2:
            # 2줄로 분할
            mid = len(words) // 2
            line1 = ' '.join(words[:mid])
            line2 = ' '.join(words[mid:])
            return f"{line1}\n{line2}"
        return label.replace('_', '\n')
    
    @staticmethod
    def _infer_node_type(label: str) -> str:
        """라벨을 기반으로 노드 타입 추론"""
        label_lower = label.lower()
        
        if any(keyword in label_lower for keyword in ['input', 'load', 'read', 'data']):
            return 'input_processing'
        elif any(keyword in label_lower for keyword in ['output', 'export', 'write', 'save']):
            return 'output'
        elif any(keyword in label_lower for keyword in ['analysis', 'count', 'calculate', 'compute']):
            return 'analysis'
        elif any(keyword in label_lower for keyword in ['grn', 'network', 'reconstruction']):
            return 'network_analysis'
        elif any(keyword in label_lower for keyword in ['trim', 'filter', 'clean', 'preprocess']):
            return 'preprocessing'
        elif any(keyword in label_lower for keyword in ['plot', 'visuali', 'graph', 'chart']):
            return 'visualization'
        else:
            return 'process'
    
    @staticmethod
    def _create_description(label: str) -> str:
        """라벨을 기반으로 설명 생성"""
        # 일반적인 설명 템플릿
        formatted_label = DOTParser._format_label(label)
        return f"Execute {formatted_label} processing step"

class D3DAGParser:
    """D3.js 호환 JSON 형식의 DAG 데이터를 파싱하는 클래스"""
    
    @staticmethod
    def parse_d3dag_json(json_content: str) -> Dict[str, Any]:
        """
        D3DAG JSON 형식을 파싱하여 노드와 엣지 정보 추출
        
        Args:
            json_content: JSON 형식의 DAG 내용
            
        Returns:
            노드와 엣지 정보를 포함한 딕셔너리
        """
        try:
            data = json.loads(json_content)
            nodes = []
            edges = []
            
            # JSON 구조 분석 및 노드/엣지 추출
            if isinstance(data, dict):
                # 노드 정보 추출
                if 'nodes' in data:
                    for node_data in data['nodes']:
                        nodes.append(D3DAGParser._process_node(node_data))
                
                # 엣지 정보 추출
                if 'links' in data or 'edges' in data:
                    edge_data = data.get('links', data.get('edges', []))
                    for edge_info in edge_data:
                        edges.append(D3DAGParser._process_edge(edge_info))
            
            return {
                'nodes': nodes,
                'edges': edges,
                'format': 'D3DAG'
            }
            
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON format: {e}")
            raise SnakemakeNativeError(
                "Failed to parse D3DAG JSON format",
                details=f"JSON decode error: {e}"
            )
        except Exception as e:
            logger.error(f"Error parsing D3DAG JSON: {e}")
            raise SnakemakeNativeError(
                "Failed to parse D3DAG format",
                details=str(e)
            )
    
    @staticmethod
    def _process_node(node_data: Dict[str, Any]) -> Dict[str, Any]:
        """노드 데이터를 처리하여 표준 형식으로 변환"""
        node_id = str(node_data.get('id', node_data.get('name', 'unknown')))
        label = node_data.get('label', node_data.get('name', node_id))
        
        return {
            'id': node_id,
            'label': DOTParser._format_label(label),
            'shortLabel': DOTParser._create_short_label(label),
            'type': DOTParser._infer_node_type(label),
            'description': DOTParser._create_description(label),
            'inputs': node_data.get('inputs', []),
            'outputs': node_data.get('outputs', []),
            'params': node_data.get('params', [])
        }
    
    @staticmethod
    def _process_edge(edge_data: Dict[str, Any]) -> Dict[str, Any]:
        """엣지 데이터를 처리하여 표준 형식으로 변환"""
        source = str(edge_data.get('source', edge_data.get('from', 'unknown')))
        target = str(edge_data.get('target', edge_data.get('to', 'unknown')))
        
        return {
            'source': source,
            'target': target,
            'label': edge_data.get('label', f'{source} → {target}')
        }

class SnakemakeNativeDAGParser:
    """
    Snakemake 네이티브 기능을 활용한 DAG 파서
    
    기존 복잡한 정규표현식 파싱을 대체하여:
    - 정확한 DAG 구조 생성
    - 유지보수성 향상
    - 성능 최적화
    """
    
    def __init__(self, timeout: int = 30):
        self.command_runner = SnakemakeCommandRunner(timeout=timeout)
    
    def parse_snakefile_with_logs(self, file_path: str, 
                                method: str = 'auto') -> Dict[str, Any]:
        """
        Snakemake 네이티브 기능을 사용하여 DAG 구조 추출 (실시간 파싱)
        
        Args:
            file_path: Snakefile 경로
            method: 파싱 방법 ('dag', 'rulegraph', 'filegraph', 'd3dag', 'auto')
            
        Returns:
            DAG 구조 정보를 포함한 딕셔너리
        """
        try:
            logger.info(f"🏗️ NATIVE PARSER: Starting real-time parse for {file_path}")
            start_time = time.time()
            
            # 파일 존재 여부 확인
            if not os.path.exists(file_path):
                logger.error(f"🏗️ NATIVE PARSER: File not found: {file_path}")
                raise SnakemakeNativeError(
                    "Snakefile not found",
                    details=f"File does not exist: {file_path}"
                )
            
            logger.error(f"🏗️ NATIVE PARSER: File exists, proceeding with parsing")
            
            # 최적 방법 결정
            if method == 'auto':
                method = self._determine_best_method(file_path)
            
            # Snakemake 명령 실행
            logger.error(f"🏗️ NATIVE PARSER: Executing Snakemake parsing with method: {method}")
            dag_data = self._execute_snakemake_parsing(file_path, method)
            logger.error(f"🏗️ NATIVE PARSER: Snakemake parsing completed, got {len(dag_data.get('nodes', []))} nodes")
            
            # 실행 순서 구축
            execution_sequence = self._build_execution_sequence(dag_data['nodes'], dag_data['edges'])
            
            # 로그 경로 정보 추가 (기존 호환성을 위해)
            logger.error(f"🏗️ NATIVE PARSER: About to call _add_log_paths for {len(dag_data['nodes'])} nodes")
            self._add_log_paths(dag_data['nodes'], file_path)
            logger.error(f"🏗️ NATIVE PARSER: _add_log_paths completed")
            
            parse_time = time.time() - start_time
            
            result = {
                'nodes': dag_data['nodes'],
                'edges': dag_data['edges'],
                'execution_sequence': execution_sequence,
                'parse_time': parse_time,
                'method': method
            }
            
            logger.info(f"Successfully parsed {len(dag_data['nodes'])} rules using {method} method in {parse_time:.3f}s")
            
            return result
            
        except SnakemakeNativeError:
            raise
        except Exception as e:
            logger.error(f"Unexpected error parsing Snakefile {file_path}: {e}")
            raise SnakemakeNativeError(
                "Unexpected parsing error",
                details=str(e)
            )
    
    def _determine_best_method(self, file_path: str) -> str:
        """파일 특성에 따라 최적의 파싱 방법 결정"""
        try:
            # 파일 크기 확인
            file_size = os.path.getsize(file_path)
            
            # 작은 파일이면 정확한 dag 방법 사용
            if file_size < 10 * 1024:  # 10KB 미만
                return 'dag'
            
            # 중간 크기 파일은 rulegraph 사용 (더 단순)
            elif file_size < 100 * 1024:  # 100KB 미만
                return 'rulegraph'
            
            # 큰 파일은 d3dag 사용 (JSON 파싱 효율성)
            else:
                return 'd3dag'
                
        except Exception:
            # 오류 시 기본 방법 사용
            return 'rulegraph'
    
    def _find_all_targets(self, file_path: str) -> List[str]:
        """Snakefile에서 모든 output 파일을 찾아서 타겟으로 반환"""
        targets = []
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            # 모든 rule의 output 찾기
            import re
            
            # rule 블록 찾기
            rules = re.split(r'^rule\s+\w+:', content, flags=re.MULTILINE)[1:]
            
            for rule_content in rules:
                # output 섹션 찾기
                output_match = re.search(r'output:\s*\n((?:\s+\w+\s*=\s*"[^"]+"\s*,?\s*\n?)+)', rule_content)
                if output_match:
                    output_section = output_match.group(1)
                    # 개별 output 파일 추출
                    file_matches = re.findall(r'\w+\s*=\s*"([^"]+)"', output_section)
                    targets.extend(file_matches)
            
            # 중복 제거하고 마지막 규칙의 output들만 사용 (일반적으로 최종 결과물)
            # 또는 모든 output 사용
            seen = set()
            unique_targets = []
            for target in targets:
                if target not in seen:
                    seen.add(target)
                    unique_targets.append(target)
                    
            # 마지막 몇 개의 타겟만 사용 (일반적으로 최종 결과물)
            if len(unique_targets) > 4:
                # 마지막 4개만 사용 (Counting_Outdegree의 4개 output)
                return unique_targets[-4:]
            
            return unique_targets
                    
        except Exception as e:
            logger.warning(f"Failed to find targets in Snakefile: {e}")
        
        return targets
    
    def _execute_snakemake_parsing(self, file_path: str, method: str) -> Dict[str, Any]:
        """Snakemake 명령을 실행하여 DAG 데이터 파싱"""
        try:
            # 타겟 찾기
            targets = self._find_all_targets(file_path)
            additional_args = targets if targets else None
            
            # 명령 실행
            stdout, stderr, returncode = self.command_runner.run_snakemake_command(
                snakefile_path=file_path,
                command_type=method,
                additional_args=additional_args
            )
            
            # 오류 확인
            if returncode != 0:
                error_msg = f"Snakemake command failed with return code {returncode}"
                logger.error(f"{error_msg}: {stderr}")
                raise SnakemakeNativeError(
                    error_msg,
                    command=f"snakemake --{method}",
                    stderr=stderr,
                    returncode=returncode
                )
            
            # 빈 출력 확인
            if not stdout.strip():
                raise SnakemakeNativeError(
                    f"Snakemake {method} command produced no output",
                    details="This may indicate an empty workflow or configuration issues"
                )
            
            # 파싱 방법에 따른 처리
            if method in ['dag', 'rulegraph', 'filegraph']:
                return DOTParser.parse_dot_graph(stdout)
            elif method == 'd3dag':
                return D3DAGParser.parse_d3dag_json(stdout)
            else:
                raise SnakemakeNativeError(f"Unsupported parsing method: {method}")
                
        except SnakemakeNativeError:
            raise
        except Exception as e:
            logger.error(f"Error executing Snakemake parsing: {e}")
            raise SnakemakeNativeError(
                "Failed to execute Snakemake parsing",
                details=str(e)
            )
    
    def _build_execution_sequence(self, nodes: List[Dict[str, Any]], 
                                edges: List[Dict[str, str]]) -> List[str]:
        """토폴로지컬 정렬을 통한 실행 순서 구축"""
        try:
            # 인접 리스트와 진입 차수 계산
            node_ids = [node['id'] for node in nodes]
            adj_list = {node_id: [] for node_id in node_ids}
            in_degree = {node_id: 0 for node_id in node_ids}
            
            for edge in edges:
                source = edge['source']
                target = edge['target']
                
                if source in adj_list and target in in_degree:
                    adj_list[source].append(target)
                    in_degree[target] += 1
            
            # 위상 정렬 (Kahn's algorithm)
            queue = [node_id for node_id, degree in in_degree.items() if degree == 0]
            sequence = []
            
            while queue:
                current = queue.pop(0)
                sequence.append(current)
                
                for neighbor in adj_list[current]:
                    in_degree[neighbor] -= 1
                    if in_degree[neighbor] == 0:
                        queue.append(neighbor)
            
            return sequence
            
        except Exception as e:
            logger.error(f"Error building execution sequence: {e}")
            # 오류 시 노드 ID 순서 반환
            return [node['id'] for node in nodes]
    
    def _add_log_paths(self, nodes: List[Dict[str, Any]], snakefile_path: str):
        """기존 API 호환성을 위해 로그 경로 정보 및 실제 input/output 파일 추가"""
        try:
            logger.error(f"🔧 PARSER DEBUG: Starting _add_log_paths for {len(nodes)} nodes")
            logger.error(f"🔧 PARSER DEBUG: Snakefile path: {snakefile_path}")
            
            # Snakefile에서 실제 rule 정보 추출
            rule_details = self._extract_rule_details_from_snakefile(snakefile_path)
            logger.error(f"🔧 PARSER DEBUG: Extracted {len(rule_details)} rule details")
            logger.error(f"🔧 PARSER DEBUG: Rule details keys: {list(rule_details.keys())}")
            
            # Snakefile 경로에서 베이스 디렉토리 추출
            base_dir = os.path.dirname(snakefile_path)
            logs_dir = os.path.join(base_dir, 'logs')
            
            for node in nodes:
                rule_id = node['id']
                rule_label = node.get('label', '').lower().replace(' ', '_')
                
                logger.error(f"🔧 PARSER DEBUG: Processing node {rule_id} with label '{node.get('label', '')}' -> '{rule_label}'")
                
                # 로그 경로 추가
                node['log_paths'] = {
                    'stdout': os.path.join(logs_dir, f'{rule_id}.stdout'),
                    'stderr': os.path.join(logs_dir, f'{rule_id}.stderr')
                }
                
                # Snakefile에서 실제 rule 정보 찾기
                rule_info = None
                
                # rule_id로 먼저 찾기 (숫자 ID인 경우)
                if rule_id in rule_details:
                    rule_info = rule_details[rule_id]
                    logger.error(f"🔧 PARSER DEBUG: Found rule_info by ID {rule_id}")
                else:
                    # rule name으로 찾기
                    for rule_name, details in rule_details.items():
                        if (rule_name.lower().replace('_', ' ') == rule_label or 
                            rule_name.lower() == rule_label or
                            details.get('label', '').lower() == rule_label):
                            rule_info = details
                            logger.error(f"🔧 PARSER DEBUG: Found rule_info by name match: {rule_name}")
                            break
                
                # rule 정보가 있으면 실제 파일 정보 추가
                if rule_info:
                    node['inputs'] = rule_info.get('inputs', [])
                    node['outputs'] = rule_info.get('outputs', [])
                    node['params'] = rule_info.get('params', [])
                    node['script'] = rule_info.get('script', "")
                    
                    logger.error(f"🔧 PARSER DEBUG: Node {rule_id} updated - inputs: {len(node['inputs'])}, outputs: {len(node['outputs'])}")
                    logger.error(f"🔧 PARSER DEBUG: Node {rule_id} outputs: {node['outputs']}")
                    
                    # 실제 로그 경로로 업데이트 (Snakefile에 정의된 경우)
                    if rule_info.get('log_paths'):
                        node['log_paths'].update(rule_info['log_paths'])
                else:
                    logger.error(f"🔧 PARSER DEBUG: No rule_info found for node {rule_id} (label: {rule_label})")
                    # 기본값 설정
                    if 'inputs' not in node:
                        node['inputs'] = []
                    if 'outputs' not in node:
                        node['outputs'] = []
                    if 'params' not in node:
                        node['params'] = []
                    if 'script' not in node:
                        node['script'] = ""
                    
        except Exception as e:
            logger.error(f"Error adding log paths and rule details: {e}")
            import traceback
            logger.error(f"Full traceback: {traceback.format_exc()}")
    
    def _extract_rule_details_from_snakefile(self, snakefile_path: str) -> Dict[str, Dict[str, Any]]:
        """Snakefile에서 실제 rule 세부정보를 추출"""
        rule_details = {}
        
        try:
            logger.error(f"🔧 EXTRACT DEBUG: Starting extraction from {snakefile_path}")
            
            # 파일 존재 확인
            if not os.path.exists(snakefile_path):
                logger.error(f"🔧 EXTRACT DEBUG: Snakefile does not exist: {snakefile_path}")
                return {}
                
            with open(snakefile_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            logger.error(f"🔧 EXTRACT DEBUG: Read {len(content)} characters from Snakefile")
            
            # rule 블록을 정규표현식으로 분리
            rule_pattern = re.compile(r'^rule\s+(\w+):(.*?)(?=^rule\s+\w+:|$)', re.MULTILINE | re.DOTALL)
            rule_matches = rule_pattern.findall(content)
            
            logger.error(f"🔧 EXTRACT DEBUG: Found {len(rule_matches)} rule matches")
            
            for i, (rule_name, rule_content) in enumerate(rule_matches):
                logger.error(f"🔧 EXTRACT DEBUG: Processing rule {i}: {rule_name}")
                
                rule_info = {
                    'name': rule_name,
                    'label': rule_name,
                    'inputs': [],
                    'outputs': [],
                    'params': [],
                    'log_paths': {},
                    'script': ""
                }
                
                # input 섹션 파싱
                input_match = re.search(r'input:\s*\n((?:\s+.*?\n)+)', rule_content)
                if input_match:
                    input_section = input_match.group(1)
                    # 파일 경로 추출 (따옴표 안의 내용)
                    file_matches = re.findall(r'["\']([^"\']+)["\']', input_section)
                    rule_info['inputs'] = [f for f in file_matches if '/' in f or '.' in f]
                    logger.error(f"🔧 EXTRACT DEBUG: Rule {rule_name} inputs: {rule_info['inputs']}")
                
                # output 섹션 파싱
                output_match = re.search(r'output:\s*\n((?:\s+.*?\n)+)', rule_content)
                if output_match:
                    output_section = output_match.group(1)
                    # 파일 경로 추출 (따옴표 안의 내용)
                    file_matches = re.findall(r'["\']([^"\']+)["\']', output_section)
                    rule_info['outputs'] = [f for f in file_matches if '/' in f or '.' in f]
                    logger.error(f"🔧 EXTRACT DEBUG: Rule {rule_name} outputs: {rule_info['outputs']}")
                else:
                    logger.error(f"🔧 EXTRACT DEBUG: No output section found for rule {rule_name}")
                
                # params 섹션 파싱
                params_match = re.search(r'params:\s*\n((?:\s+.*?\n)+)', rule_content)
                if params_match:
                    params_section = params_match.group(1)
                    # 파라미터 라인들 추출
                    param_lines = re.findall(r'\s+(\w+\s*=.*?)(?:\n|$)', params_section)
                    rule_info['params'] = param_lines
                
                # log 섹션 파싱
                log_match = re.search(r'log:\s*\n((?:\s+.*?\n)+)', rule_content)
                if log_match:
                    log_section = log_match.group(1)
                    # stdout와 stderr 로그 경로 추출
                    stdout_match = re.search(r'stdout\s*=\s*["\']([^"\']+)["\']', log_section)
                    stderr_match = re.search(r'stderr\s*=\s*["\']([^"\']+)["\']', log_section)
                    
                    if stdout_match:
                        rule_info['log_paths']['stdout'] = stdout_match.group(1)
                    if stderr_match:
                        rule_info['log_paths']['stderr'] = stderr_match.group(1)
                
                # shell 섹션 파싱 (script로 사용)
                shell_match = re.search(r'shell:\s*\n\s*["\']([^"\']+)["\']', rule_content)
                if shell_match:
                    rule_info['script'] = shell_match.group(1)
                
                # rule_name을 키로 저장
                rule_details[rule_name] = rule_info
                
                # 순차적 인덱스도 추가 (0, 1, 2, ...)
                rule_index = len(rule_details) - 1
                rule_details[str(rule_index)] = rule_info
                
                logger.error(f"🔧 EXTRACT DEBUG: Stored rule {rule_name} as index {rule_index}")
                
            logger.error(f"🔧 EXTRACT DEBUG: Final rule_details keys: {list(rule_details.keys())}")
            logger.info(f"Extracted details for {len(rule_matches)} rules from Snakefile")
            return rule_details
            
        except Exception as e:
            logger.error(f"Error extracting rule details from Snakefile: {e}")
            import traceback
            logger.error(f"Full extraction traceback: {traceback.format_exc()}")
            return {}
    
    # Cache-related methods removed for simplification

# 전역 파서 인스턴스
_native_parser = SnakemakeNativeDAGParser()

def parse_snakefile_native(file_path: str, method: str = 'auto') -> Dict[str, Any]:
    """
    편의 함수: Snakemake 네이티브 기능으로 DAG 파싱 (실시간)
    
    Args:
        file_path: Snakefile 경로
        method: 파싱 방법 ('dag', 'rulegraph', 'filegraph', 'd3dag', 'auto')
        
    Returns:
        DAG 구조 정보
    """
    return _native_parser.parse_snakefile_with_logs(file_path, method)