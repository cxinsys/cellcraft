"""진행률/타이밍/병목 집계 리포트 (PR-10 split from ``dag_parser``).

``SnakemakeRuleStatusTracker``에서 rule 상태 전체를 집계해 진행률·예상 완료 시간·병목
분석을 계산하는 메서드들을 믹스인으로 분리한다. 동작·시그니처는 원본과 동일하며, 이
믹스인은 ``SnakemakeRuleStatusTracker``에만 합성되어 사용된다.
"""
import os
import time
import logging
from typing import Dict, List, Any

# 로거 설정
logger = logging.getLogger(__name__)


class _ProgressReportMixin:
    """rule 상태 집계 기반 진행률/병목 리포트 메서드 모음 (SnakemakeRuleStatusTracker 전용 믹스인)."""

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
