"""Snakemake rule 상태 추적기 (PR-10 split from ``dag_parser``).

``SnakemakeRuleStatusTracker`` 본체 — 오케스트레이션(``get_rule_statuses``),
run.log 기반 워크플로우 완료 판정, 로그 파일 존재/내용 분석, 상태 후처리를 담는다.
개별 rule 판정 로직은 :mod:`status_analysis`, 진행률/병목 리포트는
:mod:`status_progress` 믹스인으로 분리했다. 클래스 시그니처·동작은 원본과 동일하며,
파사드가 그대로 re-export 한다.
"""
import os
import re
import time
import logging
from typing import Dict, List, Any, Optional, Tuple

from app.workflow.compiler.cache import _log_cache
from app.workflow.compiler.status_analysis import _RuleStatusAnalysisMixin
from app.workflow.compiler.status_progress import _ProgressReportMixin

# 로거 설정
logger = logging.getLogger(__name__)


class SnakemakeRuleStatusTracker(_RuleStatusAnalysisMixin, _ProgressReportMixin):
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

    def _check_logs_exist_direct(self, rule: Dict[str, Any], task_id: Optional[str] = None) -> bool:
        """
        로그 파일을 절대 경로로 직접 확인

        Args:
            rule: 룰 정보 딕셔너리
            task_id: 특정 task의 아카이브된 로그를 확인할 때 사용 (선택적)

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

            # 여러 패턴으로 로그 파일 검색
            patterns = [
                f"{rule_id}.stdout",
                f"{rule_id}.stderr",
                f"{rule_id}.log"
            ]

            # task_id가 제공되면 아카이브된 로그 먼저 확인
            if task_id:
                archived_logs_dir = os.path.join(workflow_dir, 'executions', task_id, 'logs')
                if os.path.exists(archived_logs_dir):
                    for pattern in patterns:
                        log_path = os.path.join(archived_logs_dir, pattern)
                        if os.path.exists(log_path) and os.path.getsize(log_path) > 0:
                            logger.debug(f"Found archived log file for rule {rule_id}: {log_path}")
                            return True

            # 현재 로그 디렉토리 확인 (fallback)
            logs_dir = os.path.join(workflow_dir, 'logs')

            if not os.path.exists(logs_dir):
                return False

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

    def clear_status_cache(self):
        """상태 캐시 초기화"""
        self._status_cache.clear()
        logger.debug("Rule status cache cleared")
