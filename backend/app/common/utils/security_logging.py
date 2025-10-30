"""
Security event logging module for audit trails and anomaly detection.

This module provides comprehensive security logging capabilities to track
and monitor security-relevant events across the application. It helps with:
- OWASP A09:2021 - Security Logging and Monitoring Failures
- Incident response and forensics
- Anomaly detection and pattern analysis

Events are logged in structured JSON format for easy parsing and analysis.
"""
import json
import logging
from datetime import datetime
from typing import Any, Dict, Optional
from pathlib import Path


# Configure security logger
def get_security_logger() -> logging.Logger:
    """
    Get or create the security logger instance.

    Returns a dedicated logger for security events with JSON formatting
    and appropriate handlers. The logger is configured on first access.

    Returns:
        logging.Logger: Configured security logger instance
    """
    logger = logging.getLogger("security")

    # Only configure if not already configured
    if not logger.handlers:
        # Import config lazily to avoid circular imports
        from app.common.config import SECURITY_LOG_FILE, SECURITY_LOG_LEVEL

        # Create logs directory if it doesn't exist
        log_path = Path(SECURITY_LOG_FILE)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        # Configure logger
        logger.setLevel(getattr(logging, SECURITY_LOG_LEVEL, logging.WARNING))

        # File handler for security events
        file_handler = logging.FileHandler(SECURITY_LOG_FILE, encoding='utf-8')
        file_handler.setLevel(logging.DEBUG)

        # Custom formatter that doesn't double-encode JSON
        class JSONFormatter(logging.Formatter):
            def format(self, record):
                # The message is already JSON string, just wrap it
                return json.dumps({
                    "timestamp": self.formatTime(record),
                    "level": record.levelname,
                    "message": record.getMessage()
                })

        formatter = JSONFormatter()
        file_handler.setFormatter(formatter)

        logger.addHandler(file_handler)

        # Also log CRITICAL events to console
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.CRITICAL)
        console_handler.setFormatter(logging.Formatter(
            '[SECURITY] %(levelname)s - %(message)s'
        ))
        logger.addHandler(console_handler)

        # Prevent propagation to root logger
        logger.propagate = False

    return logger


def log_security_event(
    event_type: str,
    user_id: Optional[int],
    details: Dict[str, Any],
    severity: str = "WARNING",
    ip_address: Optional[str] = None
) -> None:
    """
    Log a security event with structured data.

    Records security-relevant events including failed authentication,
    path traversal attempts, invalid file uploads, and other suspicious
    activities. All events are logged in JSON format for analysis.

    Args:
        event_type: Type of security event (e.g., "path_traversal_attempt")
        user_id: User ID associated with the event (None for unauthenticated)
        details: Dictionary containing event-specific details
        severity: Log severity level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        ip_address: Optional IP address of the requester

    Raises:
        ValueError: If severity is not a valid logging level

    Example:
        >>> log_security_event(
        ...     event_type="path_traversal_attempt",
        ...     user_id=123,
        ...     details={"attempted_path": "../../../etc/passwd"},
        ...     severity="CRITICAL"
        ... )

    Event Types:
        - path_traversal_attempt: Directory traversal attack detected
        - file_signature_mismatch: File content doesn't match extension
        - invalid_file_extension: Disallowed file type uploaded
        - file_size_exceeded: File exceeds size limit
        - empty_file_upload: Empty file uploaded
        - filename_too_long: Filename exceeds length limit
        - authentication_failure: Failed login attempt
        - authorization_failure: Access denied to resource
        - rate_limit_exceeded: Too many requests from user/IP
        - suspicious_activity: Pattern-based anomaly detected
    """
    logger = get_security_logger()

    # Validate severity
    valid_severities = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
    if severity not in valid_severities:
        raise ValueError(
            f"Invalid severity '{severity}'. Must be one of: {', '.join(valid_severities)}"
        )

    # Build structured log entry
    log_entry = {
        "event_type": event_type,
        "user_id": user_id,
        "ip_address": ip_address,
        "details": details,
        "timestamp_iso": datetime.utcnow().isoformat()
    }

    # Log at appropriate level
    log_message = json.dumps(log_entry, default=str)

    if severity == "DEBUG":
        logger.debug(log_message)
    elif severity == "INFO":
        logger.info(log_message)
    elif severity == "WARNING":
        logger.warning(log_message)
    elif severity == "ERROR":
        logger.error(log_message)
    elif severity == "CRITICAL":
        logger.critical(log_message)


def get_security_events(
    event_type: Optional[str] = None,
    user_id: Optional[int] = None,
    severity: Optional[str] = None,
    limit: int = 100
) -> list:
    """
    Retrieve recent security events for analysis.

    Reads and filters security log entries for incident investigation,
    anomaly detection, or compliance reporting. Returns events in
    reverse chronological order (most recent first).

    Args:
        event_type: Filter by specific event type (None for all)
        user_id: Filter by specific user ID (None for all)
        severity: Filter by severity level (None for all)
        limit: Maximum number of events to return (default 100)

    Returns:
        List of security event dictionaries matching filters

    Example:
        >>> events = get_security_events(
        ...     event_type="path_traversal_attempt",
        ...     severity="CRITICAL",
        ...     limit=50
        ... )
        >>> print(f"Found {len(events)} critical path traversal attempts")
    """
    from app.common.config import SECURITY_LOG_FILE

    log_path = Path(SECURITY_LOG_FILE)
    if not log_path.exists():
        return []

    events = []

    try:
        with open(SECURITY_LOG_FILE, 'r', encoding='utf-8') as f:
            for line in f:
                try:
                    # Parse JSON log entry
                    log_data = json.loads(line)

                    # Extract message field which contains the event data
                    if 'message' in log_data:
                        event = json.loads(log_data['message'])

                        # Apply filters
                        if event_type and event.get('event_type') != event_type:
                            continue
                        if user_id and event.get('user_id') != user_id:
                            continue
                        if severity and log_data.get('level') != severity:
                            continue

                        # Add log level to event
                        event['log_level'] = log_data.get('level')
                        event['log_timestamp'] = log_data.get('timestamp')

                        events.append(event)

                except json.JSONDecodeError:
                    # Skip malformed log entries
                    continue

        # Return most recent events first
        events.reverse()
        return events[:limit]

    except IOError:
        # Log file read error
        return []


def analyze_security_patterns(user_id: Optional[int] = None, time_window_hours: int = 24) -> Dict[str, Any]:
    """
    Analyze security events for patterns and anomalies.

    Performs basic anomaly detection by analyzing event patterns over a
    time window. This helps identify potential attacks, compromised accounts,
    or suspicious behavior requiring investigation.

    Args:
        user_id: Analyze patterns for specific user (None for all users)
        time_window_hours: Time window for analysis (default 24 hours)

    Returns:
        Dictionary with analysis results including:
        - total_events: Total security events in window
        - events_by_type: Count by event type
        - events_by_severity: Count by severity level
        - top_users: Users with most security events
        - anomalies: List of detected anomalies

    Example:
        >>> analysis = analyze_security_patterns(user_id=123, time_window_hours=1)
        >>> if analysis['anomalies']:
        ...     print(f"Alert: {len(analysis['anomalies'])} anomalies detected")
    """
    from datetime import timedelta

    # Get recent events
    events = get_security_events(user_id=user_id, limit=10000)

    # Filter by time window
    cutoff_time = datetime.utcnow() - timedelta(hours=time_window_hours)
    recent_events = [
        e for e in events
        if datetime.fromisoformat(e.get('timestamp_iso', '1970-01-01'))
        > cutoff_time
    ]

    # Count events by type
    events_by_type = {}
    events_by_severity = {}
    user_event_counts = {}

    for event in recent_events:
        # Count by type
        event_type = event.get('event_type', 'unknown')
        events_by_type[event_type] = events_by_type.get(event_type, 0) + 1

        # Count by severity
        severity = event.get('log_level', 'UNKNOWN')
        events_by_severity[severity] = events_by_severity.get(severity, 0) + 1

        # Count by user
        uid = event.get('user_id')
        if uid:
            user_event_counts[uid] = user_event_counts.get(uid, 0) + 1

    # Detect anomalies
    anomalies = []

    # Anomaly 1: Excessive path traversal attempts
    if events_by_type.get('path_traversal_attempt', 0) > 5:
        anomalies.append({
            "type": "excessive_path_traversal",
            "count": events_by_type['path_traversal_attempt'],
            "description": "High number of path traversal attempts detected"
        })

    # Anomaly 2: Multiple file signature mismatches
    if events_by_type.get('file_signature_mismatch', 0) > 3:
        anomalies.append({
            "type": "file_spoofing_attempts",
            "count": events_by_type['file_signature_mismatch'],
            "description": "Multiple file spoofing attempts detected"
        })

    # Anomaly 3: High rate of security events for single user
    for uid, count in user_event_counts.items():
        if count > 10:
            anomalies.append({
                "type": "user_security_events",
                "user_id": uid,
                "count": count,
                "description": f"User {uid} has {count} security events"
            })

    # Anomaly 4: High ratio of CRITICAL events
    total = len(recent_events)
    critical_count = events_by_severity.get('CRITICAL', 0)
    if total > 0 and (critical_count / total) > 0.3:
        anomalies.append({
            "type": "high_critical_ratio",
            "critical_count": critical_count,
            "total": total,
            "ratio": critical_count / total,
            "description": "High ratio of CRITICAL security events"
        })

    return {
        "total_events": len(recent_events),
        "time_window_hours": time_window_hours,
        "events_by_type": events_by_type,
        "events_by_severity": events_by_severity,
        "top_users": sorted(
            user_event_counts.items(),
            key=lambda x: x[1],
            reverse=True
        )[:10],
        "anomalies": anomalies,
        "analysis_timestamp": datetime.utcnow().isoformat()
    }
