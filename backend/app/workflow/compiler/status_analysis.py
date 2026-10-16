"""개별 rule 상태 판정 로직 (PR-10 split from ``dag_parser``).

``SnakemakeRuleStatusTracker`` (약 1,455줄)를 800줄 이하로 나누기 위해, "이 rule 하나가
pending/running/success/failed 중 무엇인가"를 결정하는 메서드들을 믹스인으로 분리한다.
동작·시그니처는 원본과 동일하며, 이 믹스인은 ``SnakemakeRuleStatusTracker``에만 합성되어
사용된다 (self 속성/다른 메서드는 합성된 클래스가 제공).
"""
import os
import logging
from typing import Dict, Any, Optional

# 로거 설정
logger = logging.getLogger(__name__)


class _RuleStatusAnalysisMixin:
    """개별 rule의 상태를 판정하는 메서드 모음 (SnakemakeRuleStatusTracker 전용 믹스인)."""

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
