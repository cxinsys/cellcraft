"""Snakefile 파싱 코어 (PR-10 split from ``dag_parser``).

``DAGParsingError`` 예외와 ``SnakemakeDAGParser`` 클래스를 담는다. 이 클래스는
Snakefile 토크나이징/rule 파싱, input-output 기반 DAG 엣지 추론, 위상 정렬 기반
실행 순서(프론트용 그래프) 직렬화를 하나의 응집된 단위로 수행한다. 기존
``dag_parser`` 모듈의 정의와 동일하며, 파사드가 그대로 re-export 한다.
"""
import os
import re
import time
import logging
from typing import Dict, List, Any, Optional
from dataclasses import dataclass

from app.workflow.compiler.cache import _dag_cache

# 로거 설정
logger = logging.getLogger(__name__)


@dataclass
class DAGParsingError(Exception):
    """DAG 파싱 관련 커스텀 예외 클래스"""
    message: str
    file_path: Optional[str] = None
    details: Optional[str] = None


class SnakemakeDAGParser:
    """
    Snakefile을 파싱하여 DAG 구조와 로그 경로를 추출하는 클래스
    """

    def __init__(self):
        self.rules = []
        self.edges = []
        self.execution_sequence = []

    def parse_snakefile_with_logs(self, file_path: str, use_cache: bool = True) -> Dict[str, Any]:
        """
        Snakefile을 파싱하여 Rule 정보와 로그 경로를 추출

        Args:
            file_path: Snakefile 경로
            use_cache: 캐시 사용 여부 (기본값: True)

        Returns:
            Dict containing nodes, edges, and execution_sequence

        Raises:
            DAGParsingError: 파일 읽기 또는 파싱 실패 시
        """
        try:
            # 캐시 확인
            if use_cache:
                cached_data = _dag_cache.get(file_path)
                if cached_data is not None:
                    logger.debug(f"Using cached DAG data for {file_path}")
                    cached_data['cached'] = True
                    return cached_data

            start_time = time.time()
            # 파일 존재 여부 확인
            if not os.path.exists(file_path):
                raise DAGParsingError(
                    f"Snakefile not found",
                    file_path=file_path,
                    details="The specified Snakefile path does not exist"
                )

            # 파일 읽기 권한 확인
            if not os.access(file_path, os.R_OK):
                raise DAGParsingError(
                    f"Cannot read Snakefile",
                    file_path=file_path,
                    details="Permission denied or file is not readable"
                )

            # 파일 크기 확인 (대용량 파일 처리)
            file_size = 0
            try:
                file_size = os.path.getsize(file_path)
                if file_size > 50 * 1024 * 1024:  # 50MB 제한
                    logger.warning(f"Large Snakefile detected: {file_size / 1024 / 1024:.1f}MB")

                if file_size > 100 * 1024 * 1024:  # 100MB 제한
                    raise DAGParsingError(
                        f"Snakefile too large",
                        file_path=file_path,
                        details=f"File size ({file_size / 1024 / 1024:.1f}MB) exceeds maximum limit (100MB)"
                    )
            except Exception as e:
                logger.warning(f"Cannot get file size for {file_path}: {e}")

            # 파일 내용 읽기 (메모리 효율적인 방법)
            try:
                content = self._read_file_efficiently(file_path, file_size)
            except UnicodeDecodeError:
                # UTF-8로 읽기 실패 시 다른 인코딩 시도
                try:
                    with open(file_path, 'r', encoding='latin-1') as f:
                        content = f.read()
                    logger.warning(f"Snakefile {file_path} read with latin-1 encoding")
                except Exception as e:
                    raise DAGParsingError(
                        f"Failed to read Snakefile with any encoding",
                        file_path=file_path,
                        details=str(e)
                    )
            except Exception as e:
                raise DAGParsingError(
                    f"Failed to read Snakefile",
                    file_path=file_path,
                    details=str(e)
                )

            # 빈 파일 확인
            if not content.strip():
                raise DAGParsingError(
                    f"Snakefile is empty",
                    file_path=file_path,
                    details="The Snakefile contains no content"
                )

            # Rule들 파싱
            rules = self.parse_rules(content)

            # Rule이 없는 경우 처리
            if not rules:
                logger.warning(f"No rules found in Snakefile: {file_path}")
                return {
                    'nodes': [],
                    'edges': [],
                    'execution_sequence': []
                }

            # 각 Rule의 로그 경로 추출
            for rule in rules:
                try:
                    rule['log_paths'] = self.extract_log_section(rule['content'])
                except Exception as e:
                    logger.error(f"Failed to extract logs for rule {rule['id']}: {e}")
                    rule['log_paths'] = {}  # 빈 로그 경로로 설정

            # 엣지 추론
            try:
                edges = self.infer_edges(rules)
            except Exception as e:
                logger.error(f"Failed to infer edges: {e}")
                edges = []  # 빈 엣지 리스트로 설정

            # 실행 순서 구축
            try:
                execution_sequence = self.build_execution_sequence(rules, edges)
            except Exception as e:
                logger.error(f"Failed to build execution sequence: {e}")
                # 규칙 ID 순서대로 기본 실행 순서 생성
                execution_sequence = [rule['id'] for rule in rules]

            parse_time = time.time() - start_time
            logger.info(f"Successfully parsed {len(rules)} rules from {file_path} in {parse_time:.3f}s")

            result = {
                'nodes': rules,
                'edges': edges,
                'execution_sequence': execution_sequence,
                'parse_time': parse_time,
                'cached': False
            }

            # 캐시에 저장
            if use_cache:
                _dag_cache.set(file_path, result)

            return result

        except DAGParsingError:
            raise
        except Exception as e:
            logger.error(f"Unexpected error parsing Snakefile {file_path}: {e}")
            raise DAGParsingError(
                f"Unexpected parsing error",
                file_path=file_path,
                details=str(e)
            )

    def parse_rules(self, content: str) -> List[Dict[str, Any]]:
        """
        Snakefile 내용에서 Rule들을 추출

        Args:
            content: Snakefile 내용

        Returns:
            List of rule dictionaries

        Raises:
            DAGParsingError: Rule 파싱 실패 시
        """
        rules = []

        try:
            # Rule 매칭 패턴: rule name: ... (다음 rule까지 또는 파일 끝까지)
            rule_regex = r'rule\s+(\w+):\s*(.*?)(?=rule\s+\w+:|$)'

            matches = re.finditer(rule_regex, content, re.DOTALL)

            for match_idx, match in enumerate(matches):
                try:
                    rule_name = match.group(1)
                    rule_content = match.group(2).strip()

                    # 빈 rule 내용 확인
                    if not rule_content:
                        logger.warning(f"Empty rule content for rule '{rule_name}'")
                        continue

                    rule_data = {
                        'id': rule_name,
                        'label': self.format_label(rule_name),
                        'shortLabel': self.create_short_label(rule_name),
                        'type': self.infer_type(rule_name),
                        'inputs': self.safe_extract_inputs(rule_content),
                        'outputs': self.safe_extract_outputs(rule_content),
                        'params': self.safe_extract_params(rule_content),
                        'description': self.create_description(rule_name),
                        'script': self.safe_extract_script(rule_content),
                        'content': rule_content  # 로그 추출을 위해 원본 내용 보존
                    }

                    rules.append(rule_data)
                    logger.debug(f"Parsed rule: {rule_name}")

                except Exception as e:
                    logger.error(f"Failed to parse rule at match {match_idx}: {e}")
                    continue

        except re.error as e:
            raise DAGParsingError(
                f"Regular expression error in rule parsing",
                details=str(e)
            )
        except Exception as e:
            raise DAGParsingError(
                f"Unexpected error parsing rules",
                details=str(e)
            )

        return rules

    def extract_log_section(self, rule_content: str) -> Dict[str, str]:
        """
        Rule 내용에서 log 섹션을 추출

        예시:
        log:
            stdout="user/testuser/workflow_1/algorithm_16/logs/TENET_Input.stdout",
            stderr="user/testuser/workflow_1/algorithm_16/logs/TENET_Input.stderr"
        """
        log_paths = {}

        # log 섹션 매칭 (다음 섹션이 나올 때까지)
        log_match = re.search(r'log:\s*(.*?)(?=\n\s*\w+:|$)', rule_content, re.DOTALL)

        if log_match:
            log_content = log_match.group(1)

            # stdout 추출
            stdout_match = re.search(r'stdout\s*=\s*["\']([^"\']+)["\']', log_content)
            if stdout_match:
                log_paths['stdout'] = stdout_match.group(1)

            # stderr 추출
            stderr_match = re.search(r'stderr\s*=\s*["\']([^"\']+)["\']', log_content)
            if stderr_match:
                log_paths['stderr'] = stderr_match.group(1)

        return log_paths

    def safe_extract_inputs(self, rule_content: str) -> List[str]:
        """Rule에서 input 파일들을 안전하게 추출"""
        try:
            return self.extract_inputs(rule_content)
        except Exception as e:
            logger.error(f"Failed to extract inputs: {e}")
            return []

    def extract_inputs(self, rule_content: str) -> List[str]:
        """Rule에서 input 파일들을 추출"""
        inputs = []

        # input 섹션 매칭
        input_match = re.search(r'input:\s*(.*?)(?=\n\s*\w+:|$)', rule_content, re.DOTALL)

        if input_match:
            input_content = input_match.group(1)

            # 파일 경로 추출 (변수명="경로" 형태)
            file_matches = re.findall(r'\w+\s*=\s*["\']([^"\']+)["\']', input_content)
            inputs.extend(file_matches)

        return inputs

    def safe_extract_outputs(self, rule_content: str) -> List[str]:
        """Rule에서 output 파일들을 안전하게 추출"""
        try:
            return self.extract_outputs(rule_content)
        except Exception as e:
            logger.error(f"Failed to extract outputs: {e}")
            return []

    def extract_outputs(self, rule_content: str) -> List[str]:
        """Rule에서 output 파일들을 추출"""
        outputs = []

        # output 섹션 매칭
        output_match = re.search(r'output:\s*(.*?)(?=\n\s*\w+:|$)', rule_content, re.DOTALL)

        if output_match:
            output_content = output_match.group(1)

            # 파일 경로 추출
            file_matches = re.findall(r'\w+\s*=\s*["\']([^"\']+)["\']', output_content)
            outputs.extend(file_matches)

        return outputs

    def safe_extract_params(self, rule_content: str) -> List[str]:
        """Rule에서 parameters를 안전하게 추출"""
        try:
            return self.extract_params(rule_content)
        except Exception as e:
            logger.error(f"Failed to extract params: {e}")
            return []

    def extract_params(self, rule_content: str) -> List[str]:
        """Rule에서 parameters를 추출"""
        params = []

        # params 섹션 매칭
        params_match = re.search(r'params:\s*(.*?)(?=\n\s*\w+:|$)', rule_content, re.DOTALL)

        if params_match:
            params_content = params_match.group(1)

            # 파라미터 추출
            param_matches = re.findall(r'(\w+)\s*=', params_content)
            params.extend(param_matches)

        return params

    def safe_extract_script(self, rule_content: str) -> str:
        """Rule에서 실행 스크립트를 안전하게 추출"""
        try:
            return self.extract_script(rule_content)
        except Exception as e:
            logger.error(f"Failed to extract script: {e}")
            return ""

    def extract_script(self, rule_content: str) -> str:
        """Rule에서 실행 스크립트 추출"""
        # shell 명령어 추출
        shell_match = re.search(r'shell:\s*["\']([^"\']+)["\']', rule_content)
        if shell_match:
            return shell_match.group(1)

        # script 파일 추출
        script_match = re.search(r'script:\s*["\']([^"\']+)["\']', rule_content)
        if script_match:
            return script_match.group(1)

        return ""

    def format_label(self, rule_name: str) -> str:
        """Rule 이름을 사용자 친화적으로 포맷팅"""
        # 언더스코어를 공백으로 변환 (원본 대소문자 유지)
        formatted = rule_name.replace('_', ' ')
        return formatted

    def create_short_label(self, rule_name: str) -> str:
        """사각형 노드 내부에 표시할 짧은 라벨 생성"""
        # 긴 이름을 줄바꿈으로 분할
        parts = rule_name.split('_')
        if len(parts) > 1:
            # 2개씩 묶어서 줄바꿈
            lines = []
            for i in range(0, len(parts), 2):
                line_parts = parts[i:i+2]
                lines.append('_'.join(line_parts))
            return '\n'.join(lines)

        return rule_name

    def infer_type(self, rule_name: str) -> str:
        """Rule 이름을 기반으로 노드 타입 추론"""
        name_lower = rule_name.lower()

        if any(keyword in name_lower for keyword in ['input', 'load', 'read']):
            return 'input_processing'
        elif any(keyword in name_lower for keyword in ['output', 'export', 'write', 'save']):
            return 'output'
        elif any(keyword in name_lower for keyword in ['analysis', 'count', 'calculate', 'compute']):
            return 'analysis'
        elif any(keyword in name_lower for keyword in ['grn', 'network', 'reconstruction']):
            return 'network_analysis'
        elif any(keyword in name_lower for keyword in ['trim', 'filter', 'clean']):
            return 'preprocessing'
        else:
            return 'process'

    def create_description(self, rule_name: str) -> str:
        """Rule 설명 생성"""
        descriptions = {
            'TENET_Input': 'Prepare input data for TENET algorithm',
            'Run': 'Execute main TENET analysis',
            'GRN_Reconstruction__FDR': 'Reconstruct gene regulatory network using FDR threshold',
            'GRN_Reconstruction__NumLinks': 'Reconstruct gene regulatory network using number of links',
            'Indirect_Edge_Trimming': 'Remove indirect edges from the network',
            'Counting_Outdegree': 'Calculate outdegree statistics for network nodes'
        }

        return descriptions.get(rule_name, f"Execute {self.format_label(rule_name)} step")

    def infer_edges(self, rules: List[Dict[str, Any]]) -> List[Dict[str, str]]:
        """Input-Output 관계를 기반으로 엣지 추론"""
        edges = []
        output_to_rule = {}

        # Output 파일과 Rule 매핑 구축
        for rule in rules:
            for output in rule['outputs']:
                normalized_output = self.normalize_filename(output)
                output_to_rule[normalized_output] = rule['id']

        # Input-Output 매칭으로 엣지 생성
        for rule in rules:
            for input_file in rule['inputs']:
                normalized_input = self.normalize_filename(input_file)

                # 매칭되는 output을 가진 rule 찾기
                for output_file, source_rule_id in output_to_rule.items():
                    if self.files_match(normalized_input, output_file):
                        edges.append({
                            'source': source_rule_id,
                            'target': rule['id'],
                            'label': os.path.basename(input_file)
                        })
                        break

        return edges

    def normalize_filename(self, filename: str) -> str:
        """파일명 정규화 (경로, 확장자 제거)"""
        # 경로에서 파일명만 추출
        basename = os.path.basename(filename)

        # 확장자 제거
        name_without_ext = os.path.splitext(basename)[0]

        return name_without_ext.lower()

    def files_match(self, input_file: str, output_file: str) -> bool:
        """입력 파일과 출력 파일이 매칭되는지 확인"""
        # 정확히 일치하는 경우
        if input_file == output_file:
            return True

        # 부분 일치 (한 파일명이 다른 파일명에 포함되는 경우)
        if input_file in output_file or output_file in input_file:
            return True

        return False

    def build_execution_sequence(self, rules: List[Dict[str, Any]], edges: List[Dict[str, str]]) -> List[str]:
        """토폴로지컬 정렬을 통한 실행 순서 구축"""
        # 인접 리스트와 진입 차수 계산
        adj_list = {rule['id']: [] for rule in rules}
        in_degree = {rule['id']: 0 for rule in rules}

        for edge in edges:
            adj_list[edge['source']].append(edge['target'])
            in_degree[edge['target']] += 1

        # 위상 정렬 (Kahn's algorithm)
        queue = [rule_id for rule_id, degree in in_degree.items() if degree == 0]
        sequence = []

        while queue:
            current = queue.pop(0)
            sequence.append(current)

            for neighbor in adj_list[current]:
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)

        return sequence

    def _read_file_efficiently(self, file_path: str, file_size: int) -> str:
        """파일을 메모리 효율적으로 읽기"""
        try:
            # 작은 파일은 한 번에 읽기
            if file_size < 10 * 1024 * 1024:  # 10MB 미만
                with open(file_path, 'r', encoding='utf-8') as f:
                    return f.read()

            # 큰 파일은 청크 단위로 읽기
            content_chunks = []
            chunk_size = 8192  # 8KB 청크

            with open(file_path, 'r', encoding='utf-8') as f:
                while True:
                    chunk = f.read(chunk_size)
                    if not chunk:
                        break
                    content_chunks.append(chunk)

                    # 메모리 사용량 모니터링 (간단한 제한)
                    if len(content_chunks) > 10000:  # 약 80MB 제한
                        logger.warning(f"File {file_path} too large for efficient processing")
                        break

            return ''.join(content_chunks)

        except MemoryError:
            raise DAGParsingError(
                f"Insufficient memory to read large Snakefile",
                file_path=file_path,
                details="File is too large to fit in memory"
            )
        except Exception as e:
            raise DAGParsingError(
                f"Error reading file efficiently",
                file_path=file_path,
                details=str(e)
            )
