"""
Unit tests for security logging and anomaly detection.

Tests the security logging system's ability to capture, store,
and analyze security events for compliance and threat detection.
"""
import json
import os
import pytest
import tempfile
from pathlib import Path
from datetime import datetime, timedelta

from app.core.logging import (
    log_security_event,
    get_security_events,
    analyze_security_patterns
)


@pytest.fixture
def temp_security_log(monkeypatch):
    """Create a temporary security log file for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        log_file = os.path.join(tmpdir, "test_security.log")

        # Monkeypatch the config to use temporary log file
        monkeypatch.setattr("app.core.config.SECURITY_LOG_FILE", log_file)

        # Reset logger handlers to force reconfiguration
        import logging
        logger = logging.getLogger("security")
        logger.handlers.clear()

        yield log_file

        # Cleanup
        logger.handlers.clear()


class TestLogSecurityEvent:
    """Test security event logging functionality."""

    def test_log_critical_event(self, temp_security_log):
        """Test logging a CRITICAL security event."""
        log_security_event(
            event_type="path_traversal_attempt",
            user_id=123,
            details={"attempted_path": "../../../etc/passwd"},
            severity="CRITICAL"
        )

        # Verify log file was created
        assert os.path.exists(temp_security_log)

        # Verify log content
        with open(temp_security_log, 'r') as f:
            log_line = f.readline()
            log_data = json.loads(log_line)

            assert log_data["level"] == "CRITICAL"
            assert "path_traversal_attempt" in log_data["message"]

            # Parse the message field
            message_data = json.loads(log_data["message"])
            assert message_data["event_type"] == "path_traversal_attempt"
            assert message_data["user_id"] == 123
            assert message_data["details"]["attempted_path"] == "../../../etc/passwd"

    def test_log_warning_event(self, temp_security_log):
        """Test logging a WARNING security event."""
        log_security_event(
            event_type="invalid_file_extension",
            user_id=456,
            details={"filename": "test.exe", "extension": ".exe"},
            severity="WARNING"
        )

        with open(temp_security_log, 'r') as f:
            log_line = f.readline()
            log_data = json.loads(log_line)
            assert log_data["level"] == "WARNING"

    def test_log_with_ip_address(self, temp_security_log):
        """Test logging with IP address."""
        log_security_event(
            event_type="authentication_failure",
            user_id=None,  # Unauthenticated
            details={"username": "admin"},
            severity="WARNING",
            ip_address="192.168.1.100"
        )

        with open(temp_security_log, 'r') as f:
            log_line = f.readline()
            log_data = json.loads(log_line)
            message_data = json.loads(log_data["message"])
            assert message_data["ip_address"] == "192.168.1.100"
            assert message_data["user_id"] is None

    def test_log_with_complex_details(self, temp_security_log):
        """Test logging with nested details."""
        log_security_event(
            event_type="file_signature_mismatch",
            user_id=789,
            details={
                "filename": "test.h5ad",
                "declared_extension": ".h5ad",
                "actual_signature": "MZ",
                "metadata": {
                    "size": 1024,
                    "uploaded_at": datetime.utcnow().isoformat()
                }
            },
            severity="CRITICAL"
        )

        with open(temp_security_log, 'r') as f:
            log_line = f.readline()
            log_data = json.loads(log_line)
            message_data = json.loads(log_data["message"])
            assert message_data["details"]["filename"] == "test.h5ad"
            assert "metadata" in message_data["details"]
            assert "size" in message_data["details"]["metadata"]

    def test_invalid_severity_raises_error(self, temp_security_log):
        """Test that invalid severity level raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            log_security_event(
                event_type="test",
                user_id=1,
                details={},
                severity="INVALID"
            )
        assert "Invalid severity" in str(exc_info.value)

    def test_multiple_events_logged(self, temp_security_log):
        """Test logging multiple events."""
        for i in range(5):
            log_security_event(
                event_type="rate_limit_exceeded",
                user_id=i,
                details={"request_count": i * 10},
                severity="WARNING"
            )

        # Verify all events were logged
        with open(temp_security_log, 'r') as f:
            lines = f.readlines()
            assert len(lines) == 5


class TestGetSecurityEvents:
    """Test security event retrieval and filtering."""

    def test_get_all_events(self, temp_security_log):
        """Test retrieving all security events."""
        # Log some events
        for i in range(3):
            log_security_event(
                event_type="test_event",
                user_id=i,
                details={"index": i},
                severity="WARNING"
            )

        # Retrieve events
        events = get_security_events()
        assert len(events) == 3
        assert all(e["event_type"] == "test_event" for e in events)

    def test_filter_by_event_type(self, temp_security_log):
        """Test filtering events by type."""
        # Log different event types
        log_security_event("path_traversal_attempt", 1, {}, "CRITICAL")
        log_security_event("invalid_file_extension", 2, {}, "WARNING")
        log_security_event("path_traversal_attempt", 3, {}, "CRITICAL")

        # Filter by type
        events = get_security_events(event_type="path_traversal_attempt")
        assert len(events) == 2
        assert all(e["event_type"] == "path_traversal_attempt" for e in events)

    def test_filter_by_user_id(self, temp_security_log):
        """Test filtering events by user ID."""
        # Log events for different users
        log_security_event("test_event", 123, {}, "WARNING")
        log_security_event("test_event", 456, {}, "WARNING")
        log_security_event("test_event", 123, {}, "WARNING")

        # Filter by user
        events = get_security_events(user_id=123)
        assert len(events) == 2
        assert all(e["user_id"] == 123 for e in events)

    def test_filter_by_severity(self, temp_security_log):
        """Test filtering events by severity."""
        # Log events with different severities
        log_security_event("test_event", 1, {}, "WARNING")
        log_security_event("test_event", 2, {}, "CRITICAL")
        log_security_event("test_event", 3, {}, "CRITICAL")

        # Filter by severity
        events = get_security_events(severity="CRITICAL")
        assert len(events) == 2
        assert all(e["log_level"] == "CRITICAL" for e in events)

    def test_limit_results(self, temp_security_log):
        """Test limiting number of results."""
        # Log many events
        for i in range(10):
            log_security_event("test_event", i, {}, "WARNING")

        # Limit to 5 events
        events = get_security_events(limit=5)
        assert len(events) == 5

    def test_events_returned_in_reverse_order(self, temp_security_log):
        """Test that events are returned most recent first."""
        # Log events with sequential IDs
        for i in range(5):
            log_security_event("test_event", i, {"order": i}, "WARNING")

        # Get events
        events = get_security_events()

        # Most recent (highest order) should be first
        assert events[0]["details"]["order"] == 4
        assert events[-1]["details"]["order"] == 0

    def test_empty_log_returns_empty_list(self, temp_security_log):
        """Test that empty log file returns empty list."""
        events = get_security_events()
        assert events == []

    def test_nonexistent_log_returns_empty_list(self, temp_security_log, monkeypatch):
        """Test that nonexistent log file returns empty list."""
        monkeypatch.setattr("app.core.config.SECURITY_LOG_FILE", "/nonexistent/path/log")
        events = get_security_events()
        assert events == []


class TestAnalyzeSecurityPatterns:
    """Test security pattern analysis and anomaly detection."""

    def test_basic_analysis(self, temp_security_log):
        """Test basic pattern analysis."""
        # Log some events
        log_security_event("path_traversal_attempt", 1, {}, "CRITICAL")
        log_security_event("invalid_file_extension", 2, {}, "WARNING")
        log_security_event("path_traversal_attempt", 1, {}, "CRITICAL")

        # Analyze patterns
        analysis = analyze_security_patterns(time_window_hours=1)

        assert analysis["total_events"] == 3
        assert "path_traversal_attempt" in analysis["events_by_type"]
        assert analysis["events_by_type"]["path_traversal_attempt"] == 2
        assert "CRITICAL" in analysis["events_by_severity"]
        assert analysis["events_by_severity"]["CRITICAL"] == 2

    def test_detect_excessive_path_traversal(self, temp_security_log):
        """Test detection of excessive path traversal attempts."""
        # Log many path traversal attempts
        for i in range(10):
            log_security_event("path_traversal_attempt", 1, {}, "CRITICAL")

        # Analyze patterns
        analysis = analyze_security_patterns()

        # Should detect anomaly
        assert len(analysis["anomalies"]) > 0
        anomaly = next(
            a for a in analysis["anomalies"]
            if a["type"] == "excessive_path_traversal"
        )
        assert anomaly["count"] == 10

    def test_detect_file_spoofing_attempts(self, temp_security_log):
        """Test detection of file spoofing attempts."""
        # Log multiple file signature mismatches
        for i in range(5):
            log_security_event("file_signature_mismatch", i, {}, "CRITICAL")

        # Analyze patterns
        analysis = analyze_security_patterns()

        # Should detect anomaly
        anomaly = next(
            (a for a in analysis["anomalies"]
             if a["type"] == "file_spoofing_attempts"),
            None
        )
        assert anomaly is not None
        assert anomaly["count"] == 5

    def test_detect_user_with_many_events(self, temp_security_log):
        """Test detection of users with excessive security events."""
        # Log many events for one user
        for i in range(15):
            log_security_event("various_events", 123, {}, "WARNING")

        # Analyze patterns
        analysis = analyze_security_patterns()

        # Should detect user anomaly
        anomaly = next(
            (a for a in analysis["anomalies"]
             if a["type"] == "user_security_events"),
            None
        )
        assert anomaly is not None
        assert anomaly["user_id"] == 123
        assert anomaly["count"] == 15

    def test_detect_high_critical_ratio(self, temp_security_log):
        """Test detection of high ratio of CRITICAL events."""
        # Log mostly CRITICAL events
        for i in range(10):
            log_security_event("critical_event", i, {}, "CRITICAL")

        # Log a few non-critical
        for i in range(2):
            log_security_event("warning_event", i, {}, "WARNING")

        # Analyze patterns
        analysis = analyze_security_patterns()

        # Should detect high critical ratio (10/12 > 0.3)
        anomaly = next(
            (a for a in analysis["anomalies"]
             if a["type"] == "high_critical_ratio"),
            None
        )
        assert anomaly is not None
        assert anomaly["critical_count"] == 10
        assert anomaly["total"] == 12

    def test_filter_by_user_in_analysis(self, temp_security_log):
        """Test analyzing patterns for specific user."""
        # Log events for different users
        for i in range(5):
            log_security_event("test_event", 123, {}, "WARNING")
        for i in range(3):
            log_security_event("test_event", 456, {}, "WARNING")

        # Analyze for specific user
        analysis = analyze_security_patterns(user_id=123)

        # Should only count events for user 123
        assert analysis["total_events"] == 5

    def test_time_window_filtering(self, temp_security_log):
        """Test that time window filtering works."""
        # This test would require manipulating event timestamps
        # For now, just verify the function accepts the parameter
        analysis = analyze_security_patterns(time_window_hours=24)
        assert analysis["time_window_hours"] == 24

    def test_top_users_sorted(self, temp_security_log):
        """Test that top users are sorted by event count."""
        # Log different numbers of events per user
        for i in range(5):
            log_security_event("test_event", 1, {}, "WARNING")
        for i in range(3):
            log_security_event("test_event", 2, {}, "WARNING")
        for i in range(10):
            log_security_event("test_event", 3, {}, "WARNING")

        # Analyze patterns
        analysis = analyze_security_patterns()

        # Top user should be user 3 with 10 events
        assert len(analysis["top_users"]) > 0
        assert analysis["top_users"][0][0] == 3
        assert analysis["top_users"][0][1] == 10

    def test_analysis_includes_timestamp(self, temp_security_log):
        """Test that analysis includes timestamp."""
        analysis = analyze_security_patterns()
        assert "analysis_timestamp" in analysis
        # Verify it's a valid ISO format timestamp
        datetime.fromisoformat(analysis["analysis_timestamp"])
