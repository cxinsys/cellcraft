import os
import re
import logging
import hashlib
import json
import time
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass
from functools import lru_cache
import threading

# 로거 설정
logger = logging.getLogger(__name__)

@dataclass
class DAGParsingError(Exception):
    """DAG 파싱 관련 커스텀 예외 클래스"""
    message: str
    file_path: Optional[str] = None
    details: Optional[str] = None

class DAGCache:
    """DAG 구조와 로그 상태를 캐싱하는 클래스"""
    def __init__(self, max_size: int = 100, ttl_seconds: int = 300):
        self.max_size = max_size
        self.ttl_seconds = ttl_seconds  # 5분 TTL
        self._cache: Dict[str, Dict[str, Any]] = {}
        self._access_times: Dict[str, float] = {}
        self._file_mtimes: Dict[str, float] = {}
        self._lock = threading.Lock()
    
    def _get_cache_key(self, file_path: str) -> str:
        """파일 경로와 수정 시간을 기반으로 캐시 키 생성"""
        try:
            file_stat = os.stat(file_path)
            mtime = file_stat.st_mtime
            size = file_stat.st_size
            
            # 파일 경로, 수정시간, 크기를 해시화
            key_data = f"{file_path}_{mtime}_{size}"
            return hashlib.md5(key_data.encode()).hexdigest()
        except (OSError, FileNotFoundError):
            # 파일이 없거나 접근할 수 없는 경우
            return hashlib.md5(file_path.encode()).hexdigest()
    
    def _is_cache_valid(self, key: str, file_path: str) -> bool:
        """캐시 유효성 검사"""
        try:
            current_time = time.time()
            
            # TTL 확인
            if key not in self._access_times:
                return False
                
            if current_time - self._access_times[key] > self.ttl_seconds:
                return False
            
            # 파일 수정 시간 확인
            if os.path.exists(file_path):
                current_mtime = os.path.getmtime(file_path)
                cached_mtime = self._file_mtimes.get(key)
                
                if cached_mtime is None or current_mtime > cached_mtime:
                    return False
            
            return True
            
        except Exception as e:
            logger.warning(f"Error checking cache validity for {file_path}: {e}")
            return False
    
    def get(self, file_path: str) -> Optional[Dict[str, Any]]:
        """캐시에서 DAG 데이터 조회"""
        with self._lock:
            key = self._get_cache_key(file_path)
            
            if key in self._cache and self._is_cache_valid(key, file_path):
                # 액세스 시간 업데이트
                self._access_times[key] = time.time()
                logger.debug(f"Cache hit for {file_path}")
                return self._cache[key].copy()
            
            # 캐시 미스 또는 무효한 캐시
            if key in self._cache:
                logger.debug(f"Cache invalidated for {file_path}")
                self._remove_key(key)
            
            return None
    
    def set(self, file_path: str, data: Dict[str, Any]) -> None:
        """DAG 데이터를 캐시에 저장"""
        with self._lock:
            key = self._get_cache_key(file_path)
            current_time = time.time()
            
            # 캐시 크기 제한 확인
            if len(self._cache) >= self.max_size:
                self._evict_oldest()
            
            # 파일 수정 시간 저장
            try:
                if os.path.exists(file_path):
                    self._file_mtimes[key] = os.path.getmtime(file_path)
            except Exception as e:
                logger.warning(f"Error getting mtime for {file_path}: {e}")
            
            self._cache[key] = data.copy()
            self._access_times[key] = current_time
            logger.debug(f"Cached DAG data for {file_path}")
    
    def _remove_key(self, key: str) -> None:
        """캐시에서 키 제거"""
        self._cache.pop(key, None)
        self._access_times.pop(key, None)
        self._file_mtimes.pop(key, None)
    
    def _evict_oldest(self) -> None:
        """가장 오래된 캐시 엔트리 제거"""
        if not self._access_times:
            return
            
        oldest_key = min(self._access_times.keys(), key=lambda k: self._access_times[k])
        self._remove_key(oldest_key)
        logger.debug(f"Evicted oldest cache entry: {oldest_key}")
    
    def clear(self) -> None:
        """전체 캐시 초기화"""
        with self._lock:
            self._cache.clear()
            self._access_times.clear()
            self._file_mtimes.clear()
            logger.info("DAG cache cleared")
    
    def get_stats(self) -> Dict[str, Any]:
        """캐시 통계 반환"""
        with self._lock:
            return {
                "size": len(self._cache),
                "max_size": self.max_size,
                "ttl_seconds": self.ttl_seconds
            }

# 전역 캐시 인스턴스
_dag_cache = DAGCache()


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


class LogFileCache:
    """로그 파일 상태를 캐싱하는 클래스"""
    def __init__(self, ttl_seconds: int = 60):  # 1분 TTL
        self.ttl_seconds = ttl_seconds
        self._cache: Dict[str, Tuple[float, bool]] = {}  # path -> (timestamp, exists)
        self._lock = threading.Lock()
    
    def file_exists(self, file_path: str) -> bool:
        """파일 존재 여부를 캐시를 통해 확인"""
        with self._lock:
            current_time = time.time()
            
            if file_path in self._cache:
                cached_time, cached_exists = self._cache[file_path]
                
                # TTL 확인
                if current_time - cached_time < self.ttl_seconds:
                    return cached_exists
            
            # 캐시 미스 또는 만료
            try:
                exists = os.path.exists(file_path)
                self._cache[file_path] = (current_time, exists)
                return exists
            except Exception as e:
                logger.warning(f"Error checking file existence {file_path}: {e}")
                return False
    
    def clear(self) -> None:
        """캐시 초기화"""
        with self._lock:
            self._cache.clear()

# 전역 로그 캐시 인스턴스
_log_cache = LogFileCache()

class SnakemakeRuleStatusTracker:
    """
    로그 파일 기반 Snakemake Rule 상태 추적 클래스
    
    향상된 기능:
    - 정확한 상태 판정 (pending, running, success, failed)
    - 로그 내용 분석을 통한 실패 감지
    - 타이밍 기반 진행률 계산
    - 스레드 안전 상태 추적
    """
    
    def __init__(self, workflow_path: str, task_status: str, use_log_cache: bool = True):
        self.workflow_path = workflow_path
        self.task_status = task_status
        self.rules = []
        self.use_log_cache = use_log_cache
        
        # 상태 추적 개선을 위한 추가 속성
        self._rule_start_patterns = {
            'job_started': ['Building DAG', 'Selecting jobs', 'Job', 'rule'],
            'rule_started': ['Submitted', 'Executing', 'Processing', 'Running'],
            'rule_finished': ['Finished job', 'Job finished', 'Complete']
        }
        
        self._error_patterns = {
            'critical': ['Error', 'CRITICAL', 'Fatal', 'Exception', 'Traceback'],
            'failure': ['Failed', 'FAILED', 'Aborted', 'terminated', 'killed'],
            'warning': ['Warning', 'WARN', 'deprecated']
        }
        
        # 상태 캐시 (룰별 상태와 마지막 확인 시간)
        self._status_cache: Dict[str, Tuple[str, float]] = {}
        self._cache_ttl = 30  # 30초 캐시
        
        # Check if workflow actually completed by analyzing run.log
        self._actual_workflow_status = self._check_workflow_completion()
        
    def get_rule_statuses(self) -> Dict[str, str]:
        """
        향상된 로그 분석을 통한 각 Rule의 상태 반환
        
        개선사항:
        - 로그 내용 분석을 통한 정확한 상태 판정
        - 타이밍 기반 상태 추론
        - 캐싱을 통한 성능 최적화
        - 더 세밀한 오류 감지
        
        Returns:
            Dict mapping rule_id to status ('pending', 'running', 'success', 'failed')
        """
        statuses = {}
        current_time = time.time()
        
        try:
            # 룰이 없는 경우 처리
            if not self.rules:
                logger.warning("No rules available for status tracking")
                return statuses
            
            # 전체 워크플로우 상태 기반 초기 판단
            workflow_completed = self.task_status in ['SUCCESS', 'FAILURE']
            workflow_failed = self.task_status == 'FAILURE'
            
            # ENHANCED DEBUG LOGGING FOR RULES
            logger.error(f"🔍 RULE DEBUG: Total rules to process: {len(self.rules)}")
            
            for i, rule in enumerate(self.rules):
                try:
                    rule_id = rule.get('id', f'unknown_rule_{i}')
                    
                    # DEBUG: Rule structure overview
                    logger.error(f"🔍 RULE DEBUG: Processing rule {i} with ID '{rule_id}'")
                    logger.error(f"🔍 RULE DEBUG: Rule {i} keys: {list(rule.keys())}")
                    
                    # 캐시 확인 (30초 TTL)
                    if rule_id in self._status_cache:
                        cached_status, cached_time = self._status_cache[rule_id]
                        if current_time - cached_time < self._cache_ttl:
                            statuses[rule_id] = cached_status
                            continue
                    
                    # 향상된 상태 분석
                    status = self._analyze_rule_status(rule, i, workflow_completed, workflow_failed)
                    
                    # 캐시 업데이트
                    self._status_cache[rule_id] = (status, current_time)
                    statuses[rule_id] = status
                    
                except Exception as e:
                    rule_id = rule.get('id', f'unknown_rule_{i}')
                    logger.error(f"Error processing rule {rule_id}: {e}")
                    statuses[rule_id] = 'pending'  # 기본값으로 설정
            
            # 후처리: 상태 일관성 검사 및 조정
            statuses = self._post_process_statuses(statuses, workflow_completed, workflow_failed)
            
        except Exception as e:
            logger.error(f"Unexpected error in get_rule_statuses: {e}")
            # 기본 상태로 모든 룰을 pending으로 설정
            for i, rule in enumerate(self.rules):
                rule_id = rule.get('id', f'unknown_rule_{i}')
                statuses[rule_id] = 'pending'
        
        return statuses
    
    def _check_workflow_completion(self) -> Optional[str]:
        """
        Check if the workflow has actually completed by analyzing run.log
        
        Returns:
            'SUCCESS' if workflow completed successfully (100% done)
            'FAILURE' if workflow failed
            None if status cannot be determined
        """
        try:
            # Construct path to run.log
            if not self.workflow_path:
                return None
                
            workflow_dir = os.path.dirname(self.workflow_path)
            run_log_path = os.path.join(workflow_dir, 'logs', 'run.log')
            
            if not os.path.exists(run_log_path):
                logger.debug(f"run.log not found at {run_log_path}")
                return None
            
            # Read last few lines of run.log to check completion
            with open(run_log_path, 'r') as f:
                lines = f.readlines()
                
            # Check for completion patterns
            for line in reversed(lines[-20:]):  # Check last 20 lines
                # Success patterns
                if '100% done' in line or 'Complete log:' in line:
                    logger.info(f"Workflow completed successfully based on run.log")
                    return 'SUCCESS'
                    
                # Failure patterns
                if 'Error' in line or 'Failed' in line or 'Traceback' in line:
                    # But ignore if it's part of the completion message
                    if 'Complete log:' not in line:
                        logger.info(f"Workflow failed based on run.log")
                        return 'FAILURE'
            
            # Check for "X of Y steps (100%) done" pattern
            import re
            for line in reversed(lines[-50:]):  # Check more lines for this pattern
                match = re.search(r'(\d+)\s+of\s+(\d+)\s+steps\s+\((\d+)%\)\s+done', line)
                if match:
                    completed = int(match.group(1))
                    total = int(match.group(2))
                    percentage = int(match.group(3))
                    
                    if percentage == 100 and completed == total:
                        logger.info(f"Workflow completed successfully: {completed}/{total} steps (100%) done")
                        return 'SUCCESS'
                    else:
                        logger.debug(f"Workflow progress: {completed}/{total} steps ({percentage}%) done")
                        
            return None
            
        except Exception as e:
            logger.error(f"Error checking workflow completion: {e}")
            return None
    
    def safe_check_log_existence(self, log_files: List[str]) -> Dict[str, Any]:
        """로그 파일 존재 여부와 크기를 안전하게 확인 (캐싱 지원)"""
        try:
            result = {
                "has_any_logs": False,
                "log_details": {},
                "total_files": len([f for f in log_files if f]),
                "existing_files": 0
            }
            
            if not log_files:
                return result
                
            for log_path in log_files:
                if not log_path:  # 빈 경로 무시
                    continue
                
                log_info = {
                    "exists": False,
                    "size": 0,
                    "accessible": False,
                    "error": None
                }
                
                try:
                    # 캐시 사용 여부에 따라 분기
                    if self.use_log_cache:
                        exists = _log_cache.file_exists(log_path)
                    else:
                        exists = os.path.exists(log_path)
                    
                    log_info["exists"] = exists
                    
                    if exists:
                        try:
                            # 파일 크기 확인
                            file_stat = os.stat(log_path)
                            log_info["size"] = file_stat.st_size
                            log_info["accessible"] = True
                            result["existing_files"] += 1
                            
                            # 의미있는 크기의 로그 파일이 있으면 True
                            if file_stat.st_size > 0:
                                result["has_any_logs"] = True
                            
                        except (OSError, PermissionError) as e:
                            log_info["error"] = str(e)
                            logger.warning(f"Cannot access log file stats {log_path}: {e}")
                        
                except (OSError, PermissionError) as e:
                    log_info["error"] = str(e)
                    logger.warning(f"Cannot check log file {log_path}: {e}")
                    continue
                
                # 파일명만 키로 사용 (경로가 너무 길어지지 않도록)
                log_filename = os.path.basename(log_path) if log_path else "unknown"
                result["log_details"][log_filename] = log_info
                        
            return result
            
        except Exception as e:
            logger.error(f"Error checking log files existence: {e}")
            return {
                "has_any_logs": False,
                "log_details": {},
                "total_files": 0,
                "existing_files": 0,
                "error": str(e)
            }
    
    def _analyze_rule_status(self, rule: Dict[str, Any], rule_index: int, 
                           workflow_completed: bool, workflow_failed: bool) -> str:
        """
        Robust rule status analysis with prioritized logic
        
        Priority Order:
        1. Output files existence (most reliable)
        2. Actual workflow status from run.log
        3. Task status with workflow completion context
        4. Fallback to log-based analysis
        
        Args:
            rule: 룰 정보 딕셔너리
            rule_index: 룰의 실행 순서 인덱스
            workflow_completed: 워크플로우 완료 여부 (deprecated)
            workflow_failed: 워크플로우 실패 여부 (deprecated)
            
        Returns:
            상태 문자열 ('pending', 'running', 'success', 'failed')
        """
        try:
            rule_id = rule.get('id', f'unknown_rule_{rule_index}')
            
            # ENHANCED DEBUG LOGGING
            logger.error(f"🔍 DEBUG: Analyzing rule {rule_id} (index {rule_index})")
            logger.error(f"🔍 DEBUG: Rule object keys: {list(rule.keys())}")
            logger.error(f"🔍 DEBUG: Rule outputs: {rule.get('outputs', [])}")
            logger.error(f"🔍 DEBUG: Rule log_paths: {rule.get('log_paths', {})}")
            logger.error(f"🔍 DEBUG: Workflow path: {self.workflow_path}")
            logger.error(f"🔍 DEBUG: Task status: {self.task_status}, Actual workflow status: {self._actual_workflow_status}")
            
            logger.debug(f"Analyzing rule {rule_id} (index {rule_index}) - task_status: {self.task_status}, actual_workflow_status: {self._actual_workflow_status}")
            
            # PRIORITY 1: 출력 파일 존재 확인 (가장 확실한 완료 증거)
            if self._check_rule_outputs_exist(rule):
                logger.info(f"Rule {rule_id} -> SUCCESS (output files exist)")
                return 'success'
            
            # PRIORITY 2: 실제 워크플로우 완료 상태 확인
            if self._actual_workflow_status == 'SUCCESS':
                # 워크플로우가 100% 완료된 경우, outputs나 로그 파일이 있으면 성공으로 간주
                outputs = rule.get('outputs', [])
                log_paths = rule.get('log_paths', {})
                
                if outputs:
                    # 출력 파일 정의가 있으면 이 룰이 실행되었음을 의미 
                    logger.error(f"✅ PRIORITY 2A: Rule {rule_id} -> SUCCESS (workflow completed + outputs defined: {len(outputs)} files)")
                    return 'success'
                elif log_paths:
                    # 로그 경로가 정의되어 있으면 이 룰이 워크플로우의 일부임을 의미
                    logger.error(f"✅ PRIORITY 2B: Rule {rule_id} -> SUCCESS (workflow completed + log paths defined)")
                    return 'success'
                elif self._has_any_execution_evidence(rule):
                    logger.info(f"Rule {rule_id} -> SUCCESS (workflow completed + execution evidence)")
                    return 'success'
                else:
                    # 출력도 로그 경로도 없고 실행 증거도 없으면 이 룰은 실제로 실행되지 않았을 수 있음
                    logger.debug(f"Rule {rule_id} -> PENDING (workflow completed but no evidence)")
                    return 'pending'
            
            # PRIORITY 2.5: RUNNING 상태에서 outputs가 있는 rule 처리
            if self.task_status == 'RUNNING':
                outputs = rule.get('outputs', [])
                if outputs:
                    # RUNNING 상태에서 outputs가 정의된 rule은 실행 순서를 기반으로 상태 결정
                    status = self._determine_running_rule_status(rule, rule_index)
                    logger.error(f"✅ PRIORITY 2.5: Rule {rule_id} -> {status.upper()} (RUNNING + outputs defined: {len(outputs)} files)")
                    return status
            
            # PRIORITY 3: 실제 워크플로우 실패 상태 처리
            if self._actual_workflow_status == 'FAILURE':
                return self._handle_failed_workflow_rule(rule, rule_index)
            
            # PRIORITY 4: 태스크 상태 기반 처리 (REVOKED 등)
            effective_status = self.task_status
            
            if effective_status == 'SUCCESS':
                # DB 상태가 SUCCESS이지만 실제 워크플로우 상태가 불분명한 경우
                if self._has_any_execution_evidence(rule):
                    logger.info(f"Rule {rule_id} -> SUCCESS (task success + execution evidence)")
                    return 'success'
                else:
                    logger.debug(f"Rule {rule_id} -> PENDING (task success but no execution evidence)")
                    return 'pending'
                    
            elif effective_status in ['FAILURE', 'REVOKED']:
                # 실제 워크플로우가 완료되었다면 성공으로 처리 (cancellation after completion case)
                if self._actual_workflow_status == 'SUCCESS':
                    if self._has_any_execution_evidence(rule):
                        if effective_status == 'REVOKED':
                            logger.info(f"Rule {rule_id} -> SUCCESS (workflow completed before task was revoked)")
                        else:
                            logger.info(f"Rule {rule_id} -> SUCCESS (workflow completed before task failed)")
                        return 'success'
                return self._handle_failed_task_rule(rule, rule_index)
            
            # PRIORITY 5: 실행 중인 워크플로우에 대한 기존 로직 (마지막 수단)
            return self._analyze_running_rule_status(rule, rule_index)
                        
        except Exception as e:
            logger.error(f"Error analyzing rule status for {rule.get('id', 'unknown')}: {e}")
            return 'pending'
    
    def _determine_running_rule_status(self, rule: Dict[str, Any], rule_index: int) -> str:
        """
        RUNNING 상태에서 rule의 상태를 결정
        실행 순서와 로그 파일 존재 여부를 기반으로 판단
        
        Args:
            rule: 룰 정보 딕셔너리
            rule_index: 룰 인덱스
            
        Returns:
            'success', 'running', 또는 'pending'
        """
        try:
            rule_id = rule.get('id', str(rule_index))
            
            # 1. 로그 파일 존재 확인으로 완료 여부 판단
            if self._check_logs_exist_direct(rule):
                logger.error(f"🏃 RUNNING DEBUG: Rule {rule_id} has log files - likely completed")
                return 'success'
            
            # 2. 실행 순서 기반 추론
            # 이전 rule들이 모두 완료되었고, 현재 rule에 대한 증거가 없으면 pending
            # 현재 rule이 마지막 실행된 rule이면 running
            last_running_rule_index = self._find_last_running_rule_index()
            
            if last_running_rule_index is not None:
                if rule_index < last_running_rule_index:
                    # 마지막 실행 rule보다 이전 - 완료된 것으로 추정
                    logger.error(f"🏃 RUNNING DEBUG: Rule {rule_id} (index {rule_index}) before last running ({last_running_rule_index}) - SUCCESS")
                    return 'success'
                elif rule_index == last_running_rule_index:
                    # 마지막 실행 rule - 현재 실행 중
                    logger.error(f"🏃 RUNNING DEBUG: Rule {rule_id} (index {rule_index}) is last running - RUNNING")
                    return 'running'
                else:
                    # 마지막 실행 rule 이후 - 아직 대기 중
                    logger.error(f"🏃 RUNNING DEBUG: Rule {rule_id} (index {rule_index}) after last running ({last_running_rule_index}) - PENDING")
                    return 'pending'
            
            # 3. 기본값: 첫 번째 rule은 running, 나머지는 pending
            if rule_index == 0:
                logger.error(f"🏃 RUNNING DEBUG: Rule {rule_id} is first rule - RUNNING")
                return 'running'
            else:
                logger.error(f"🏃 RUNNING DEBUG: Rule {rule_id} default - PENDING")
                return 'pending'
                
        except Exception as e:
            logger.error(f"Error determining running rule status: {e}")
            return 'pending'
    
    def _has_any_execution_evidence(self, rule: Dict[str, Any]) -> bool:
        """
        룰이 실행되었다는 증거가 있는지 확인 (다중 검증)
        
        Args:
            rule: 룰 정보 딕셔너리
            
        Returns:
            실행 증거가 있으면 True, 없으면 False
        """
        try:
            # 1. 출력 파일 존재 (가장 확실한 증거)
            if self._check_rule_outputs_exist(rule):
                return True
            
            # 2. 로그 파일 존재 (절대 경로로 직접 확인)
            if self._check_logs_exist_direct(rule):
                return True
            
            # 3. 기존 로그 경로 방식으로 확인 (fallback)
            log_paths = rule.get('log_paths', {})
            if log_paths:
                log_files = [path for path in log_paths.values() if path]
                log_check_result = self.safe_check_log_existence(log_files)
                if log_check_result.get("has_any_logs", False):
                    return True
            
            # 4. 룰 내용이 있고 의미있는 스크립트가 정의됨
            rule_content = rule.get('content', '')
            if rule_content and rule_content.strip() and len(rule_content.strip()) > 10:
                # 단순히 내용이 있다고 실행된 것은 아니지만, 다른 증거와 함께 고려
                return False  # 이것만으로는 실행 증거가 아님
                
            return False
            
        except Exception as e:
            logger.debug(f"Error checking execution evidence for rule {rule.get('id', 'unknown')}: {e}")
            return False
    
    def _check_logs_exist_direct(self, rule: Dict[str, Any]) -> bool:
        """
        로그 파일을 절대 경로로 직접 확인
        
        Args:
            rule: 룰 정보 딕셔너리
            
        Returns:
            로그 파일이 존재하면 True, 없으면 False
        """
        try:
            rule_id = rule.get('id', '')
            if not rule_id:
                return False
            
            # 워크플로우 디렉토리에서 로그 파일 직접 검색
            workflow_dir = os.path.dirname(self.workflow_path) if self.workflow_path else None
            if not workflow_dir:
                return False
                
            logs_dir = os.path.join(workflow_dir, 'logs')
            
            if not os.path.exists(logs_dir):
                return False
            
            # 여러 패턴으로 로그 파일 검색
            patterns = [
                f"{rule_id}.stdout",
                f"{rule_id}.stderr", 
                f"{rule_id}.log"
            ]
            
            for pattern in patterns:
                log_path = os.path.join(logs_dir, pattern)
                if os.path.exists(log_path) and os.path.getsize(log_path) > 0:
                    logger.debug(f"Found log file for rule {rule_id}: {log_path}")
                    return True
            
            logger.debug(f"No log files found for rule {rule_id} in {logs_dir}")
            return False
            
        except Exception as e:
            logger.debug(f"Error checking direct logs for rule {rule.get('id', 'unknown')}: {e}")
            return False
    
    def _handle_failed_workflow_rule(self, rule: Dict[str, Any], rule_index: int) -> str:
        """
        실패한 워크플로우에서 룰 상태를 처리 (태스크 실패 처리와 일관성 유지)
        
        Args:
            rule: 룰 정보 딕셔너리
            rule_index: 룰의 실행 순서 인덱스
            
        Returns:
            상태 문자열 ('pending', 'success', 'failed')
        """
        try:
            rule_id = rule.get('id', f'unknown_rule_{rule_index}')
            
            # 마지막으로 실행된 rule을 찾아서 실패 처리
            last_running_rule_index = self._find_last_running_rule_index()
            
            if last_running_rule_index is not None:
                if rule_index == last_running_rule_index:
                    # 마지막 실행 rule -> 실패한 rule
                    logger.error(f"💥 WORKFLOW FAILURE: Rule {rule_id} -> FAILED (last running rule when workflow failed)")
                    return 'failed'
                elif rule_index < last_running_rule_index:
                    # 마지막 실행 rule 이전의 rule들 중 실행 증거가 있는 것은 완료된 것으로 처리
                    if self._has_any_execution_evidence(rule):
                        logger.info(f"Rule {rule_id} -> SUCCESS (completed before workflow failure)")
                        return 'success'
                    else:
                        logger.debug(f"Rule {rule_id} -> PENDING (no execution evidence before failure)")
                        return 'pending'
                else:
                    # 마지막 실행 rule 이후의 rule들은 실행되지 않음
                    logger.debug(f"Rule {rule_id} -> PENDING (after failed rule)")
                    return 'pending'
            else:
                # 마지막 실행 rule을 찾을 수 없는 경우 - 기존 로직 사용
                if self._has_any_execution_evidence(rule):
                    logger.info(f"Rule {rule_id} -> SUCCESS (has execution evidence, no clear failure point)")
                    return 'success'
                else:
                    logger.debug(f"Rule {rule_id} -> PENDING (no execution evidence, no clear failure point)")
                    return 'pending'
                
        except Exception as e:
            logger.error(f"Error handling failed workflow rule {rule.get('id', 'unknown')}: {e}")
            return 'pending'
    
    def _handle_failed_task_rule(self, rule: Dict[str, Any], rule_index: int) -> str:
        """
        실패한 태스크에서 룰 상태를 처리 (기존 로직 유지)
        
        Args:
            rule: 룰 정보 딕셔너리
            rule_index: 룰의 실행 순서 인덱스
            
        Returns:
            상태 문자열 ('pending', 'success', 'failed')
        """
        try:
            rule_id = rule.get('id', f'unknown_rule_{rule_index}')
            
            # Find the last rule that was running when task failed
            last_running_rule_index = self._find_last_running_rule_index()
            
            if last_running_rule_index is not None:
                if rule_index == last_running_rule_index:
                    task_action = "revoked" if self.task_status == 'REVOKED' else "failed"
                    if self.task_status == 'REVOKED':
                        logger.error(f"🚫 TASK REVOKED: Rule {rule_id} -> FAILED (last running rule when task was revoked)")
                    else:
                        logger.error(f"💥 TASK FAILED: Rule {rule_id} -> FAILED (last running rule when task failed)")
                    return 'failed'
                elif rule_index < last_running_rule_index:
                    # Rules before the failed rule that have execution evidence are considered completed
                    if self._has_any_execution_evidence(rule):
                        logger.info(f"Rule {rule_id} -> SUCCESS (completed before failure)")
                        return 'success'
                    else:
                        return 'pending'
                else:
                    # Rules after the failed rule remain pending
                    logger.debug(f"Rule {rule_id} -> PENDING (after failed rule)")
                    return 'pending'
            else:
                # No clear last running rule found, use execution evidence
                if self._has_any_execution_evidence(rule):
                    logger.info(f"Rule {rule_id} -> SUCCESS (has execution evidence, no clear failure point)")
                    return 'success'
                else:
                    logger.debug(f"Rule {rule_id} -> PENDING (no execution evidence, no clear failure point)")
                    return 'pending'
                
        except Exception as e:
            logger.error(f"Error handling failed task rule {rule.get('id', 'unknown')}: {e}")
            return 'pending'
    
    def _analyze_running_rule_status(self, rule: Dict[str, Any], rule_index: int) -> str:
        """
        실행 중인 워크플로우의 룰 상태 분석 (기존 로직 간소화)
        
        Args:
            rule: 룰 정보 딕셔너리
            rule_index: 룰의 실행 순서 인덱스
            
        Returns:
            상태 문자열 ('pending', 'running', 'success')
        """
        try:
            rule_id = rule.get('id', f'unknown_rule_{rule_index}')
            
            # 실행 증거가 있는지 확인
            if self._has_any_execution_evidence(rule):
                # 로그 내용을 분석해서 완료/실행중 구분
                log_paths = rule.get('log_paths', {})
                if log_paths:
                    log_analysis = self._analyze_log_content(log_paths)
                    
                    # Check for explicit failure patterns
                    if log_analysis.get('has_critical_error', False) or log_analysis.get('has_failure', False):
                        logger.info(f"Rule {rule_id} -> FAILED (error patterns in logs)")
                        return 'failed'
                    
                    # Check for explicit completion patterns
                    if log_analysis.get('has_completion_marker', False):
                        logger.info(f"Rule {rule_id} -> SUCCESS (completion markers in logs)")
                        return 'success'
                
                # 증거는 있지만 명확한 완료/실패 패턴이 없으면 실행중으로 간주
                logger.info(f"Rule {rule_id} -> RUNNING (has execution evidence)")
                return 'running'
            else:
                # 실행 증거가 없으면 아직 시작되지 않음
                logger.debug(f"Rule {rule_id} -> PENDING (no execution evidence)")
                return 'pending'
                
        except Exception as e:
            logger.error(f"Error analyzing running rule status for {rule.get('id', 'unknown')}: {e}")
            return 'pending'

    def _find_last_running_rule_index(self) -> Optional[int]:
        """
        Find the index of the last rule that was running when the task failed
        
        Returns:
            Index of the last running rule, or None if not found
        """
        try:
            for i in range(len(self.rules) - 1, -1, -1):  # Reverse order
                rule = self.rules[i]
                log_paths = rule.get('log_paths', {})
                log_files = [path for path in log_paths.values() if path]
                log_check_result = self.safe_check_log_existence(log_files)
                
                # If this rule has logs, it was at least started
                if log_check_result.get("has_any_logs", False):
                    # Check if it completed or was still running
                    if self._check_rule_outputs_exist(rule):
                        # If outputs exist, it completed, continue looking
                        continue
                    else:
                        # No outputs but has logs - this was likely the running rule
                        return i
            
            return None
            
        except Exception as e:
            logger.error(f"Error finding last running rule: {e}")
            return None
    
    def _check_rule_outputs_exist(self, rule: Dict[str, Any]) -> bool:
        """
        강화된 룰 출력 파일 존재 확인
        
        Args:
            rule: Rule dictionary containing outputs information
            
        Returns:
            True if expected outputs exist, False otherwise
        """
        try:
            outputs = rule.get('outputs', [])
            
            # ENHANCED DEBUG LOGGING
            logger.error(f"📁 DEBUG: Checking outputs for rule {rule.get('id', 'unknown')}")
            logger.error(f"📁 DEBUG: Outputs from rule: {outputs}")
            logger.error(f"📁 DEBUG: Workflow path: {self.workflow_path}")
            
            if not outputs:
                # No outputs defined, cannot determine completion from outputs
                logger.error(f"📁 DEBUG: No outputs defined for rule {rule.get('id', 'unknown')}")
                return False
            
            workflow_dir = os.path.dirname(self.workflow_path) if self.workflow_path else None
            logger.error(f"📁 DEBUG: Workflow dir: {workflow_dir}")
            
            for output in outputs:
                if not output:
                    continue
                    
                # Handle both string and list outputs
                output_path = output
                if isinstance(output, dict):
                    output_path = output.get('path') or output.get('file')
                    if not output_path:
                        continue
                
                # Convert relative paths to absolute paths with Docker support
                if not os.path.isabs(output_path):
                    # Docker 컨테이너에서는 working directory가 다를 수 있음
                    base_dirs = [
                        '/app',  # Docker 컨테이너 내부
                        '/home/dmshin/cellcraft/backend',  # 호스트 시스템
                        os.getcwd(),  # 현재 작업 디렉토리
                    ]
                    
                    if workflow_dir:
                        base_dirs.insert(0, workflow_dir)
                    
                    possible_paths = []
                    for base_dir in base_dirs:
                        possible_paths.extend([
                            os.path.join(base_dir, output_path),
                            os.path.join(base_dir, 'backend', output_path) if not output_path.startswith('backend') else os.path.join(base_dir, output_path),
                        ])
                    
                    # Add original path as fallback
                    possible_paths.append(output_path)
                else:
                    possible_paths = [output_path]
                
                logger.error(f"📁 DEBUG: Trying paths for output {output_path}: {possible_paths[:3]}...")  # Show first 3 paths
                
                # Check if any of the possible paths exist
                for path in possible_paths:
                    try:
                        if os.path.exists(path):
                            file_size = os.path.getsize(path)
                            logger.error(f"📁 DEBUG: Found file {path} with size {file_size}")
                            if file_size > 0:
                                logger.error(f"✅ SUCCESS: Rule {rule.get('id', 'unknown')} has valid output: {path}")
                                return True
                    except (OSError, IOError) as e:
                        logger.error(f"📁 DEBUG: Error checking path {path}: {e}")
                        continue
            
            logger.debug(f"Rule {rule.get('id', 'unknown')} - no outputs found")
            return False
            
        except Exception as e:
            logger.error(f"Error checking rule outputs for {rule.get('id', 'unknown')}: {e}")
            return False
    
    def _is_subsequent_rule_running(self, current_rule_index: int) -> bool:
        """
        Check if any subsequent rule is running, indicating current rule has completed
        
        Args:
            current_rule_index: Index of current rule in execution sequence
            
        Returns:
            True if a subsequent rule is running, False otherwise
        """
        try:
            # Check rules after current rule
            for i in range(current_rule_index + 1, len(self.rules)):
                subsequent_rule = self.rules[i]
                log_paths = subsequent_rule.get('log_paths', {})
                log_files = [path for path in log_paths.values() if path]
                log_check_result = self.safe_check_log_existence(log_files)
                
                # If subsequent rule has logs, it means it started
                if log_check_result.get("has_any_logs", False):
                    # Check if current rule's outputs are inputs to this subsequent rule
                    if self._is_dependency_relationship(current_rule_index, i):
                        logger.debug(f"Rule {current_rule_index} completed because rule {i} is running")
                        return True
                    # Even without direct dependency, if next rule started, current likely completed
                    elif i == current_rule_index + 1:  # Immediate next rule
                        logger.debug(f"Rule {current_rule_index} completed because next rule {i} is running")
                        return True
            
            return False
            
        except Exception as e:
            logger.error(f"Error checking subsequent rules: {e}")
            return False
    
    def _is_dependency_relationship(self, current_rule_index: int, subsequent_rule_index: int) -> bool:
        """
        Check if there's a dependency relationship between two rules
        
        Args:
            current_rule_index: Index of the current rule
            subsequent_rule_index: Index of the subsequent rule
            
        Returns:
            True if subsequent rule depends on current rule's outputs
        """
        try:
            if current_rule_index >= len(self.rules) or subsequent_rule_index >= len(self.rules):
                return False
                
            current_rule = self.rules[current_rule_index]
            subsequent_rule = self.rules[subsequent_rule_index]
            
            current_outputs = current_rule.get('outputs', [])
            subsequent_inputs = subsequent_rule.get('inputs', [])
            
            # Check if any output of current rule is input to subsequent rule
            for output in current_outputs:
                for input_file in subsequent_inputs:
                    # Normalize paths for comparison
                    output_normalized = os.path.basename(output) if output else ""
                    input_normalized = os.path.basename(input_file) if input_file else ""
                    
                    if output_normalized and input_normalized:
                        if output_normalized == input_normalized:
                            return True
                        # Also check if one is contained in the other (partial matches)
                        if output_normalized in input_normalized or input_normalized in output_normalized:
                            return True
            
            return False
            
        except Exception as e:
            logger.error(f"Error checking dependency relationship: {e}")
            return False
    
    def _analyze_log_content(self, log_paths: Dict[str, str]) -> Dict[str, bool]:
        """
        로그 파일 내용을 분석하여 상태 판정에 필요한 정보 추출 (개선된 패턴 감지)
        
        Returns:
            분석 결과 딕셔너리:
            - has_start_marker: 시작 표시 있음
            - has_completion_marker: 완료 표시 있음
            - has_critical_error: 심각한 오류 있음
            - has_failure: 실패 표시 있음
            - has_warning: 경고 있음
        """
        analysis = {
            'has_start_marker': False,
            'has_completion_marker': False,
            'has_critical_error': False,
            'has_failure': False,
            'has_warning': False
        }
        
        try:
            # 향상된 패턴 정의
            enhanced_start_patterns = [
                'job', 'rule', 'executing', 'processing', 'running', 'started',
                'building dag', 'selecting jobs', 'submitted', 'begin', 'starting'
            ]
            
            enhanced_completion_patterns = [
                'finished job', 'job finished', 'complete', 'completed', 'done',
                'success', 'successful', 'finished', 'end of job', 'job completed'
            ]
            
            enhanced_error_patterns = [
                'error:', 'critical:', 'fatal:', 'exception:', 'traceback',
                'failed:', 'failure:', 'aborted', 'terminated', 'killed',
                'could not', 'cannot', 'unable to', 'permission denied',
                'file not found', 'no such file', 'command not found'
            ]
            
            # stdout 로그 분석 (더 포괄적)
            if 'stdout' in log_paths and log_paths['stdout']:
                stdout_content = self._read_log_file(log_paths['stdout'])
                if stdout_content:
                    content_lower = stdout_content.lower()
                    
                    # 시작 패턴 확인 (개선)
                    analysis['has_start_marker'] = any(
                        pattern in content_lower for pattern in enhanced_start_patterns
                    )
                    
                    # 완료 패턴 확인 (개선)
                    analysis['has_completion_marker'] = any(
                        pattern in content_lower for pattern in enhanced_completion_patterns
                    )
                    
                    # stdout에서도 오류 확인 (일부 도구는 stdout에 오류 출력)
                    if any(pattern in content_lower for pattern in enhanced_error_patterns):
                        analysis['has_critical_error'] = True
            
            # stderr 로그 분석 (오류 패턴)
            if 'stderr' in log_paths and log_paths['stderr']:
                stderr_content = self._read_log_file(log_paths['stderr'])
                if stderr_content:
                    content_lower = stderr_content.lower()
                    
                    # 내용이 있으면 일단 시작된 것으로 간주 (stderr에 뭔가 기록됨)
                    if content_lower.strip():
                        analysis['has_start_marker'] = True
                    
                    # 심각한 에러 패턴 확인
                    analysis['has_critical_error'] = any(
                        pattern in content_lower for pattern in enhanced_error_patterns
                    )
                    
                    # 실패 패턴은 critical error와 동일하게 처리
                    analysis['has_failure'] = analysis['has_critical_error']
                    
                    # 경고 패턴 확인
                    warning_patterns = ['warning', 'warn', 'deprecated', 'notice']
                    analysis['has_warning'] = any(
                        pattern in content_lower for pattern in warning_patterns
                    )
        
        except Exception as e:
            logger.error(f"Error analyzing log content: {e}")
        
        return analysis
    
    def _read_log_file(self, log_path: str, max_size: int = 5 * 1024 * 1024) -> Optional[str]:
        """
        로그 파일을 안전하게 읽기
        
        Args:
            log_path: 로그 파일 경로
            max_size: 최대 읽을 파일 크기 (기본 5MB)
            
        Returns:
            로그 파일 내용 또는 None
        """
        try:
            # 파일 존재 확인
            file_exists = _log_cache.file_exists(log_path) if self.use_log_cache else os.path.exists(log_path)
            if not file_exists:
                return None
            
            # 읽기 권한 확인
            if not os.access(log_path, os.R_OK):
                logger.warning(f"Cannot read log file: {log_path}")
                return None
            
            # 파일 크기 확인
            file_size = os.path.getsize(log_path)
            
            if file_size == 0:
                return ""
            
            if file_size > max_size:
                logger.warning(f"Log file too large ({file_size} bytes): {log_path}")
                # 앞부분과 끝부분 읽기
                with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
                    # 처음 1KB
                    start_content = f.read(1024)
                    # 마지막 2KB
                    f.seek(max(0, file_size - 2048))
                    end_content = f.read()
                    return start_content + "\n...\n" + end_content
            else:
                with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
                    return f.read()
            
        except Exception as e:
            logger.error(f"Error reading log file {log_path}: {e}")
            return None
    
    def _post_process_statuses(self, statuses: Dict[str, str], 
                             workflow_completed: bool, workflow_failed: bool) -> Dict[str, str]:
        """
        상태 일관성 검사 및 조정
        
        Args:
            statuses: 룰별 상태 딕셔너리
            workflow_completed: 워크플로우 완료 여부
            workflow_failed: 워크플로우 실패 여부
            
        Returns:
            조정된 상태 딕셔너리
        """
        try:
            # 1. 워크플로우 완료 시 일관성 검사
            if workflow_completed:
                if workflow_failed:
                    # 실패한 워크플로우에서 마지막 실패 룰 찾기
                    last_failed = self._find_last_failed_rule(statuses)
                    if last_failed:
                        statuses[last_failed] = 'failed'
                        # 실패한 룰 이후의 모든 룰은 pending
                        self._mark_subsequent_rules_pending(statuses, last_failed)
                else:
                    # 성공한 워크플로우에서 모든 실행된 룰은 success
                    for rule_id, status in statuses.items():
                        if status in ['running']:
                            statuses[rule_id] = 'success'
            
            # 2. RUNNING 상태 워크플로우에서 마지막 완료된 단계를 running으로 표시
            if self.task_status == 'RUNNING':
                statuses = self._ensure_last_non_pending_is_running(statuses)
            
            # 3. running 상태 검증 (완료된 워크플로우에는 running이 없어야 함)
            if workflow_completed:
                for rule_id, status in statuses.items():
                    if status == 'running':
                        if workflow_failed:
                            statuses[rule_id] = 'failed'
                        else:
                            statuses[rule_id] = 'success'
            
        except Exception as e:
            logger.error(f"Error in post-processing statuses: {e}")
        
        return statuses
    
    def _ensure_last_non_pending_is_running(self, statuses: Dict[str, str]) -> Dict[str, str]:
        """
        RUNNING 상태에서 진행된 단계들(pending이 아닌) 중 마지막을 running으로 표시
        
        Args:
            statuses: 현재 상태 딕셔너리
            
        Returns:
            수정된 상태 딕셔너리
        """
        try:
            if not self.rules:
                return statuses
                
            # 실행 순서대로 룰을 확인하여 마지막으로 완료된 룰 찾기
            last_non_pending_index = -1
            last_non_pending_rule_id = None
            
            for i, rule in enumerate(self.rules):
                rule_id = rule.get('id', f'unknown_rule_{i}')
                current_status = statuses.get(rule_id, 'pending')
                
                # pending이 아닌 상태 (success, running, failed) 중 마지막을 찾기
                if current_status != 'pending':
                    last_non_pending_index = i
                    last_non_pending_rule_id = rule_id
                    
            logger.info(f"🏃 RUNNING ENHANCEMENT: Found last non-pending rule at index {last_non_pending_index}: {last_non_pending_rule_id}")
                    
            # 마지막으로 완료된 룰을 running으로 변경
            if last_non_pending_rule_id and last_non_pending_index >= 0:
                original_status = statuses.get(last_non_pending_rule_id)
                statuses[last_non_pending_rule_id] = 'running'
                logger.info(f"🏃 RUNNING ENHANCEMENT: Changed rule {last_non_pending_rule_id} from '{original_status}' to 'running'")
                
                # 만약 다른 running 상태가 있다면 success로 변경 (한 번에 하나만 running)
                for rule_id, status in statuses.items():
                    if rule_id != last_non_pending_rule_id and status == 'running':
                        statuses[rule_id] = 'success'
                        logger.info(f"🏃 RUNNING ENHANCEMENT: Changed previous running rule {rule_id} to 'success'")
            else:
                logger.info("🏃 RUNNING ENHANCEMENT: No non-pending rules found, keeping original statuses")
                
        except Exception as e:
            logger.error(f"Error in _ensure_last_non_pending_is_running: {e}")
            
        return statuses
    
    def _find_last_failed_rule(self, statuses: Dict[str, str]) -> Optional[str]:
        """
        실패한 룰 중 마지막 룰을 찾기
        """
        try:
            # 실행 순서 역순으로 확인
            for rule in reversed(self.rules):
                rule_id = rule.get('id', 'unknown_rule')
                if rule_id in statuses and statuses[rule_id] in ['running', 'failed']:
                    # 로그를 다시 확인하여 실제로 오류가 있는지 검사
                    log_paths = rule.get('log_paths', {})
                    log_analysis = self._analyze_log_content(log_paths)
                    
                    if log_analysis['has_critical_error'] or log_analysis['has_failure']:
                        return rule_id
                    
            return None
        except Exception as e:
            logger.error(f"Error finding last failed rule: {e}")
            return None
    
    def _mark_subsequent_rules_pending(self, statuses: Dict[str, str], failed_rule_id: str):
        """
        실패한 룰 이후의 모든 룰을 pending으로 표시
        """
        try:
            failed_rule_found = False
            
            for rule in self.rules:
                rule_id = rule.get('id', 'unknown_rule')
                
                if rule_id == failed_rule_id:
                    failed_rule_found = True
                    continue
                
                if failed_rule_found and rule_id in statuses:
                    # 실패한 룰 이후의 룰들은 pending으로 설정
                    if statuses[rule_id] not in ['pending']:
                        statuses[rule_id] = 'pending'
                        
        except Exception as e:
            logger.error(f"Error marking subsequent rules pending: {e}")
    
    def has_error_in_logs(self, log_paths: Dict[str, str]) -> bool:
        """
        로그 파일에서 오류 여부를 안전하게 확인 (하위 호환성 유지)
        """
        try:
            log_analysis = self._analyze_log_content(log_paths)
            return log_analysis['has_critical_error'] or log_analysis['has_failure']
        except Exception as e:
            logger.error(f"Error checking error logs: {e}")
            return False
    
    def find_last_executed_rule(self) -> Optional[str]:
        """Failure 시 마지막으로 실행된 Rule을 안전하게 찾기 (하위 호환성 유지)"""
        try:
            return self._find_last_failed_rule({})
        except Exception as e:
            logger.error(f"Error in find_last_executed_rule: {e}")
            return None
    
    def get_enhanced_progress_info(self) -> Dict[str, Any]:
        """
        향상된 진행률 정보 계산
        
        Returns:
            상세한 진행률 정보를 포함한 딕셔너리:
            - basic_progress: 기본 진행률 정보
            - timing_info: 타이밍 기반 정보
            - estimated_completion: 예상 완료 시간
            - bottleneck_analysis: 병목 분석
        """
        try:
            statuses = self.get_rule_statuses()
            
            # 기본 진행률 계산
            basic_progress = self._calculate_basic_progress(statuses)
            
            # 타이밍 분석 (로그 파일 수정 시간 기반)
            timing_info = self._analyze_timing_info(statuses)
            
            # 예상 완료 시간 계산
            estimated_completion = self._estimate_completion_time(basic_progress, timing_info)
            
            # 병목 분석
            bottleneck_analysis = self._analyze_bottlenecks(statuses, timing_info)
            
            return {
                "basic_progress": basic_progress,
                "timing_info": timing_info,
                "estimated_completion": estimated_completion,
                "bottleneck_analysis": bottleneck_analysis,
                "rule_statuses": statuses
            }
            
        except Exception as e:
            logger.error(f"Error calculating enhanced progress info: {e}")
            return self._get_fallback_progress_info()
    
    def _calculate_basic_progress(self, statuses: Dict[str, str]) -> Dict[str, Any]:
        """기본 진행률 계산"""
        try:
            total_rules = len(statuses)
            if total_rules == 0:
                return {
                    "total_rules": 0,
                    "completed_rules": 0,
                    "failed_rules": 0,
                    "running_rules": 0,
                    "pending_rules": 0,
                    "percentage": 0.0,
                    "completion_percentage": 0.0,
                    "is_stalled": False
                }
            
            # 상태별 카운트
            completed_rules = sum(1 for status in statuses.values() if status == 'success')
            failed_rules = sum(1 for status in statuses.values() if status == 'failed')
            running_rules = sum(1 for status in statuses.values() if status == 'running')
            pending_rules = sum(1 for status in statuses.values() if status == 'pending')
            
            # 진행률 계산 (성공 + 실패를 완료로 간주)
            finished_rules = completed_rules + failed_rules
            percentage = (finished_rules / total_rules * 100) if total_rules > 0 else 0
            
            # 성공률 (성공한 것만)
            completion_percentage = (completed_rules / total_rules * 100) if total_rules > 0 else 0
            
            # 정체 상태 감지 (running이지만 너무 오래 지속)
            is_stalled = self._detect_stall_condition(statuses)
            
            return {
                "total_rules": total_rules,
                "completed_rules": completed_rules,
                "failed_rules": failed_rules,
                "running_rules": running_rules,
                "pending_rules": pending_rules,
                "percentage": round(percentage, 1),
                "completion_percentage": round(completion_percentage, 1),
                "is_stalled": is_stalled
            }
        
        except Exception as e:
            logger.error(f"Error calculating basic progress: {e}")
            return {"total_rules": 0, "percentage": 0.0, "error": str(e)}
    
    def _analyze_timing_info(self, statuses: Dict[str, str]) -> Dict[str, Any]:
        """로그 파일 수정 시간을 기반으로 타이밍 정보 분석"""
        try:
            timing_data = []
            current_time = time.time()
            
            for i, rule in enumerate(self.rules):
                try:
                    rule_id = rule.get('id', f'unknown_rule_{i}')
                    log_paths = rule.get('log_paths', {})
                    
                    if rule_id not in statuses:
                        continue
                    
                    # 로그 파일 타임스탬프 수집
                    log_timestamps = []
                    for log_type, log_path in log_paths.items():
                        if log_path and os.path.exists(log_path):
                            try:
                                mtime = os.path.getmtime(log_path)
                                log_timestamps.append({
                                    "log_type": log_type,
                                    "timestamp": mtime,
                                    "age_seconds": current_time - mtime
                                })
                            except OSError:
                                continue
                    
                    if log_timestamps:
                        # 가장 최근 로그 타임스탬프 사용
                        latest_log = max(log_timestamps, key=lambda x: x["timestamp"])
                        
                        timing_data.append({
                            "rule_id": rule_id,
                            "rule_index": i,
                            "status": statuses[rule_id],
                            "latest_log_time": latest_log["timestamp"],
                            "age_seconds": latest_log["age_seconds"],
                            "log_count": len(log_timestamps)
                        })
                        
                except Exception as e:
                    logger.warning(f"Error analyzing timing for rule {rule_id}: {e}")
                    continue
            
            # 타이밍 통계 계산
            if timing_data:
                ages = [t["age_seconds"] for t in timing_data]
                avg_age = sum(ages) / len(ages)
                
                # 실행 시간 추정 (연속된 룰들의 시간 차이)
                execution_times = []
                sorted_timing = sorted(timing_data, key=lambda x: x["rule_index"])
                
                for i in range(1, len(sorted_timing)):
                    if (sorted_timing[i]["status"] in ['success', 'failed'] and 
                        sorted_timing[i-1]["status"] in ['success', 'failed']):
                        
                        time_diff = sorted_timing[i-1]["latest_log_time"] - sorted_timing[i]["latest_log_time"]
                        if 0 < time_diff < 3600:  # 1시간 이내의 합리적인 시간 차이만
                            execution_times.append(abs(time_diff))
                
                avg_execution_time = sum(execution_times) / len(execution_times) if execution_times else 0
                
                return {
                    "total_rules_with_logs": len(timing_data),
                    "average_log_age_seconds": round(avg_age, 1),
                    "estimated_avg_rule_time": round(avg_execution_time, 1),
                    "timing_data": timing_data[:10],  # 최대 10개만 반환
                    "has_timing_data": True
                }
            else:
                return {"has_timing_data": False, "message": "No timing data available"}
                
        except Exception as e:
            logger.error(f"Error analyzing timing info: {e}")
            return {"has_timing_data": False, "error": str(e)}
    
    def _estimate_completion_time(self, basic_progress: Dict[str, Any], 
                                timing_info: Dict[str, Any]) -> Dict[str, Any]:
        """예상 완료 시간 계산"""
        try:
            if not timing_info.get("has_timing_data", False):
                return {"available": False, "message": "Insufficient timing data"}
            
            total_rules = basic_progress.get("total_rules", 0)
            completed_rules = basic_progress.get("completed_rules", 0)
            failed_rules = basic_progress.get("failed_rules", 0)
            remaining_rules = total_rules - completed_rules - failed_rules
            
            if remaining_rules <= 0:
                return {
                    "available": True,
                    "estimated_seconds_remaining": 0,
                    "estimated_completion_time": "Completed",
                    "confidence": "high"
                }
            
            # 평균 실행 시간 기반 예상
            avg_rule_time = timing_info.get("estimated_avg_rule_time", 0)
            
            if avg_rule_time > 0:
                estimated_seconds = remaining_rules * avg_rule_time
                
                # 신뢰도 계산
                confidence = "high" if len(timing_info.get("timing_data", [])) >= 3 else "low"
                
                current_time = time.time()
                completion_time = current_time + estimated_seconds
                
                return {
                    "available": True,
                    "estimated_seconds_remaining": round(estimated_seconds, 1),
                    "estimated_completion_time": time.strftime('%Y-%m-%d %H:%M:%S', 
                                                             time.localtime(completion_time)),
                    "confidence": confidence,
                    "remaining_rules": remaining_rules,
                    "avg_rule_time_seconds": avg_rule_time
                }
            else:
                return {"available": False, "message": "Cannot estimate execution time"}
                
        except Exception as e:
            logger.error(f"Error estimating completion time: {e}")
            return {"available": False, "error": str(e)}
    
    def _analyze_bottlenecks(self, statuses: Dict[str, str], 
                           timing_info: Dict[str, Any]) -> Dict[str, Any]:
        """병목 현상 분석"""
        try:
            bottlenecks = []
            
            # running 상태가 너무 오래 지속되는 룰 찾기
            if timing_info.get("has_timing_data", False):
                timing_data = timing_info.get("timing_data", [])
                avg_rule_time = timing_info.get("estimated_avg_rule_time", 60)
                
                for timing in timing_data:
                    if timing["status"] == "running":
                        # 평균 실행 시간의 2배 이상 소요되면 병목으로 간주
                        if timing["age_seconds"] > max(avg_rule_time * 2, 300):  # 최소 5분
                            bottlenecks.append({
                                "rule_id": timing["rule_id"],
                                "type": "long_running",
                                "duration_seconds": timing["age_seconds"],
                                "expected_max_seconds": avg_rule_time * 2,
                                "severity": "high" if timing["age_seconds"] > avg_rule_time * 5 else "medium"
                            })
            
            # 연속된 실패 패턴 찾기
            consecutive_failures = 0
            for rule in self.rules:
                rule_id = rule.get('id', 'unknown_rule')
                if rule_id in statuses:
                    if statuses[rule_id] == 'failed':
                        consecutive_failures += 1
                    else:
                        if consecutive_failures >= 2:
                            bottlenecks.append({
                                "type": "consecutive_failures",
                                "failure_count": consecutive_failures,
                                "severity": "critical"
                            })
                        consecutive_failures = 0
            
            return {
                "has_bottlenecks": len(bottlenecks) > 0,
                "bottleneck_count": len(bottlenecks),
                "bottlenecks": bottlenecks,
                "recommendations": self._generate_bottleneck_recommendations(bottlenecks)
            }
            
        except Exception as e:
            logger.error(f"Error analyzing bottlenecks: {e}")
            return {"has_bottlenecks": False, "error": str(e)}
    
    def _generate_bottleneck_recommendations(self, bottlenecks: List[Dict[str, Any]]) -> List[str]:
        """병목 현상에 대한 권장사항 생성"""
        recommendations = []
        
        try:
            for bottleneck in bottlenecks:
                if bottleneck["type"] == "long_running":
                    recommendations.append(
                        f"Rule '{bottleneck['rule_id']}' has been running for "
                        f"{bottleneck['duration_seconds']:.0f}s. Consider checking resource availability."
                    )
                elif bottleneck["type"] == "consecutive_failures":
                    recommendations.append(
                        f"Multiple consecutive failures detected. "
                        f"Check input data and configuration."
                    )
            
            if not recommendations:
                recommendations.append("No specific recommendations at this time.")
                
        except Exception as e:
            logger.error(f"Error generating recommendations: {e}")
            recommendations = ["Unable to generate recommendations due to error."]
        
        return recommendations
    
    def _detect_stall_condition(self, statuses: Dict[str, str]) -> bool:
        """정체 상태 감지"""
        try:
            running_count = sum(1 for status in statuses.values() if status == 'running')
            
            # running 상태인 룰이 있는지 확인
            if running_count == 0:
                return False
            
            # 타이밍 정보 기반 정체 감지 (간단한 버전)
            current_time = time.time()
            
            for rule in self.rules:
                rule_id = rule.get('id', 'unknown_rule')
                if rule_id in statuses and statuses[rule_id] == 'running':
                    log_paths = rule.get('log_paths', {})
                    
                    for log_path in log_paths.values():
                        if log_path and os.path.exists(log_path):
                            try:
                                mtime = os.path.getmtime(log_path)
                                age = current_time - mtime
                                
                                # 30분 이상 로그 업데이트가 없으면 정체로 간주
                                if age > 1800:  # 30분
                                    return True
                            except OSError:
                                continue
            
            return False
            
        except Exception as e:
            logger.error(f"Error detecting stall condition: {e}")
            return False
    
    def _get_fallback_progress_info(self) -> Dict[str, Any]:
        """오류 시 기본 진행률 정보 반환"""
        try:
            basic_statuses = {}
            for i, rule in enumerate(self.rules):
                rule_id = rule.get('id', f'unknown_rule_{i}')
                basic_statuses[rule_id] = 'pending'
            
            return {
                "basic_progress": {
                    "total_rules": len(self.rules),
                    "percentage": 0.0,
                    "error": "Could not calculate detailed progress"
                },
                "timing_info": {"has_timing_data": False},
                "estimated_completion": {"available": False},
                "bottleneck_analysis": {"has_bottlenecks": False},
                "rule_statuses": basic_statuses
            }
        except Exception as e:
            return {
                "error": f"Failed to generate fallback progress info: {e}",
                "basic_progress": {"total_rules": 0, "percentage": 0.0}
            }


    def clear_status_cache(self):
        """상태 캐시 초기화"""
        self._status_cache.clear()
        logger.debug("Rule status cache cleared")


# 캐시 관리 유틸리티 함수들
def clear_all_caches():
    """모든 캐시를 초기화"""
    _dag_cache.clear()
    _log_cache.clear()
    logger.info("All caches cleared")

def get_cache_stats() -> Dict[str, Any]:
    """캐시 통계 정보 반환"""
    return {
        "dag_cache": _dag_cache.get_stats(),
        "log_cache": {
            "size": len(_log_cache._cache),
            "ttl_seconds": _log_cache.ttl_seconds
        }
    }

def set_cache_config(dag_max_size: int = None, dag_ttl: int = None, log_ttl: int = None):
    """캐시 설정 변경"""
    global _dag_cache, _log_cache
    
    if dag_max_size is not None:
        _dag_cache.max_size = dag_max_size
        
    if dag_ttl is not None:
        _dag_cache.ttl_seconds = dag_ttl
        
    if log_ttl is not None:
        _log_cache.ttl_seconds = log_ttl
    
    logger.info(f"Cache config updated: DAG(max_size={_dag_cache.max_size}, ttl={_dag_cache.ttl_seconds}s), Log(ttl={_log_cache.ttl_seconds}s)")

@lru_cache(maxsize=128)
def get_file_hash(file_path: str, mtime: float) -> str:
    """파일 해시 계산 (LRU 캐시 적용)"""
    try:
        hash_data = f"{file_path}_{mtime}"
        return hashlib.md5(hash_data.encode()).hexdigest()
    except Exception as e:
        logger.error(f"Error calculating file hash for {file_path}: {e}")
        return hashlib.md5(file_path.encode()).hexdigest()