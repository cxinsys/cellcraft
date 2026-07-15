"""Visualization Snakefile 생성/추출 (PR-10 split from ``workflow/utils``).

compile 파이프라인의 일부인 Snakefile rule 추출·파라미터 치환·독립 실행 Snakefile
생성 로직을 ``workflow/utils.py``에서 이곳(compiler 하위)으로 이동했다. 함수
시그니처·동작은 원본과 동일하며, ``app.workflow.utils``가 그대로 re-export 하므로
기존 참조 경로(``app.workflow.utils.extract_rule_block`` 등)는 유지된다.
"""
import os
import re
import logging
from typing import Dict, Optional, Tuple, List, Any

from app.core.exceptions import SnakefileGenerationError
from app.workflow.errors import (
    FileValidationError,
    SnakefileProcessingError,
)

# Configure logging
logger = logging.getLogger(__name__)


def extract_rule_block(snakefile_path: str, visualization_title: str, result_path: str,
                      parameters: Optional[Dict[str, Any]] = None,
                      file_mappings: Optional[Dict[str, str]] = None) -> Tuple[str, str]:
    """
    Enhanced rule extraction with dynamic parameter replacement and better error handling.

    Args:
        snakefile_path: Path to the Snakefile containing rules
        visualization_title: Name of the rule to extract
        result_path: Path where extracted rule should be saved
        parameters: Optional dict of parameters to replace in the rule
        file_mappings: Optional dict mapping file placeholders to actual paths

    Returns:
        Tuple of (extracted_content, result_path)

    Raises:
        SnakefileProcessingError: If snakefile processing fails
        FileValidationError: If files don't exist or are invalid
    """
    logger.info(f"Extracting rule '{visualization_title}' from {snakefile_path}")

    # Validate inputs
    if not snakefile_path:
        raise SnakefileProcessingError(
            "Snakefile path is required",
            snakefile_path=snakefile_path,
            rule_name=visualization_title
        )

    if not visualization_title:
        raise SnakefileProcessingError(
            "Rule name is required",
            snakefile_path=snakefile_path,
            rule_name=visualization_title
        )

    if not os.path.exists(snakefile_path):
        raise FileValidationError(
            f"Snakefile not found",
            file_path=snakefile_path,
            context={"operation": "extract_rule_block", "rule": visualization_title}
        )

    if not os.path.isfile(snakefile_path):
        raise FileValidationError(
            f"Path is not a file",
            file_path=snakefile_path,
            context={"operation": "extract_rule_block", "rule": visualization_title}
        )

    # rule 블록들을 분리하기 위한 시작 패턴
    rule_pattern = f"rule {visualization_title}:"

    # 파일 내용을 줄 단위로 분리
    try:
        with open(snakefile_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
    except UnicodeDecodeError as e:
        raise FileValidationError(
            f"Failed to read Snakefile: encoding error",
            file_path=snakefile_path,
            context={"error": str(e), "encoding": "utf-8"}
        )
    except Exception as e:
        raise SnakefileProcessingError(
            f"Failed to read Snakefile: {str(e)}",
            snakefile_path=snakefile_path,
            rule_name=visualization_title
        ) from e

    result = []
    is_target_rule = False
    current_indent = 0
    rule_found = False

    for line_num, line in enumerate(lines, 1):
        # 대상 rule 블록의 시작을 찾음
        if line.startswith(rule_pattern):
            is_target_rule = True
            rule_found = True
            current_indent = len(line) - len(line.lstrip())
            result.append(line)
            logger.debug(f"Found rule '{visualization_title}' at line {line_num}")
            continue

        # rule 블록 내부의 내용을 수집
        if is_target_rule:
            # 현재 줄의 들여쓰기 수준 확인
            if line.strip() == "":
                result.append(line)  # Preserve empty lines within the rule
                continue

            line_indent = len(line) - len(line.lstrip())

            # 들여쓰기가 같거나 큰 경우에만 해당 rule에 속한 것으로 간주
            if line_indent > current_indent:
                result.append(line)
            else:
                break

    if not rule_found:
        try:
            available_rules = _get_available_rules(snakefile_path)
        except Exception:
            available_rules = ["Unable to determine available rules"]
        raise SnakefileProcessingError(
            f"Rule '{visualization_title}' not found. Available rules: {available_rules}",
            snakefile_path=snakefile_path,
            rule_name=visualization_title
        )

    # 빈 줄 제거하고 결과 문자열 생성
    result = [line.rstrip() for line in result if line.strip() or line.strip() == ""]
    extracted_content = 'import os\n' + '\n'.join(result)

    # Apply parameter and file mapping replacements if provided
    if parameters or file_mappings:
        extracted_content = _apply_replacements(extracted_content, parameters, file_mappings)

    # 추출한 내용을 파일에 쓰기
    try:
        # Ensure parent directory exists
        result_dir = os.path.dirname(result_path)
        if result_dir:
            os.makedirs(result_dir, exist_ok=True)

        with open(result_path, 'w', encoding='utf-8') as file:
            file.write(extracted_content)
        logger.info(f"Successfully extracted rule to {result_path}")
    except PermissionError as e:
        raise FileValidationError(
            f"Permission denied writing to result path",
            file_path=result_path,
            context={"operation": "write_extracted_rule", "error": str(e)}
        )
    except Exception as e:
        raise SnakefileProcessingError(
            f"Failed to write extracted rule: {str(e)}",
            snakefile_path=snakefile_path,
            rule_name=visualization_title
        ) from e

    return extracted_content, result_path


def generate_visualization_snakefile(plugin_path: str, selected_script_name: str,
                                   visualization_params: Dict[str, Any],
                                   file_mappings: Dict[str, str],
                                   output_path: str) -> str:
    """
    Generate a standalone executable Snakefile for a specific visualization script.

    This function reads a plugin's Snakefile, extracts only the rule for the selected script,
    replaces parameter placeholders with actual values, and generates a standalone Snakefile
    that can be executed independently.

    Args:
        plugin_path: Path to the plugin directory containing Snakefile
        selected_script_name: Name of the visualization rule/script to extract
        visualization_params: Dictionary of visualization parameters and their values
        file_mappings: Dictionary mapping input file placeholders to actual file paths
        output_path: Path where the generated Snakefile should be saved

    Returns:
        str: Path to the generated Snakefile

    Raises:
        SnakefileProcessingError: If Snakefile generation fails
        FileValidationError: If plugin or files don't exist

    Example:
        snakefile_path = generate_visualization_snakefile(
            plugin_path="/path/to/GRNViz",
            selected_script_name="Heatmap",
            visualization_params={"Top_Genes": 20, "Sample_Size": 150},
            file_mappings={"expression.csv": "/data/expr.csv", "input.txt": "/data/network.txt"},
            output_path="/tmp/heatmap_snakefile.smk"
        )
    """
    logger.info(f"Generating visualization Snakefile for '{selected_script_name}' from {plugin_path}")

    # Validate inputs
    if not plugin_path:
        raise SnakefileProcessingError(
            "Plugin path is required",
            rule_name=selected_script_name
        )

    if not selected_script_name:
        raise SnakefileProcessingError(
            "Script name is required",
            snakefile_path=plugin_path
        )

    if not output_path:
        raise SnakefileProcessingError(
            "Output path is required",
            snakefile_path=plugin_path,
            rule_name=selected_script_name
        )

    # Validate plugin path
    if not os.path.exists(plugin_path):
        raise FileValidationError(
            f"Plugin path not found",
            file_path=plugin_path,
            context={"script_name": selected_script_name, "operation": "generate_snakefile"}
        )

    if not os.path.isdir(plugin_path):
        raise FileValidationError(
            f"Plugin path is not a directory",
            file_path=plugin_path,
            context={"script_name": selected_script_name}
        )

    # Construct Snakefile path
    snakefile_path = os.path.join(plugin_path, "Snakefile")
    if not os.path.exists(snakefile_path):
        raise FileValidationError(
            f"Snakefile not found in plugin",
            file_path=snakefile_path,
            context={"plugin_path": plugin_path, "script_name": selected_script_name}
        )

    # Validate the rule exists before extraction
    if not validate_visualization_rule(snakefile_path, selected_script_name):
        try:
            available_rules = _get_available_rules(snakefile_path)
        except Exception:
            available_rules = ["Unable to determine available rules"]
        raise SnakefileProcessingError(
            f"Visualization rule '{selected_script_name}' not found. Available rules: {available_rules}",
            snakefile_path=snakefile_path,
            rule_name=selected_script_name
        )

    # Extract the specific rule with parameter replacements
    try:
        extracted_content, _ = extract_rule_block(
            snakefile_path=snakefile_path,
            visualization_title=selected_script_name,
            result_path=output_path,
            parameters=visualization_params,
            file_mappings=file_mappings
        )

        logger.info(f"Successfully generated visualization Snakefile: {output_path}")
        return output_path

    except (SnakefileProcessingError, FileValidationError):
        # Re-raise our custom exceptions
        raise
    except Exception as e:
        raise SnakefileProcessingError(
            f"Failed to generate visualization Snakefile: {str(e)}",
            snakefile_path=snakefile_path,
            rule_name=selected_script_name
        ) from e


def validate_visualization_rule(snakefile_path: str, rule_name: str) -> bool:
    """
    Check if a specific rule name exists in a Snakefile and validate its structure.

    Args:
        snakefile_path: Path to the Snakefile to check
        rule_name: Name of the rule to validate

    Returns:
        bool: True if rule exists and has valid structure, False otherwise

    Raises:
        FileNotFoundError: If snakefile_path doesn't exist
    """
    logger.debug(f"Validating rule '{rule_name}' in {snakefile_path}")

    if not os.path.exists(snakefile_path):
        raise FileNotFoundError(f"Snakefile not found: {snakefile_path}")

    try:
        with open(snakefile_path, 'r', encoding='utf-8') as file:
            content = file.read()

        # Check if rule exists
        rule_pattern = f"rule {rule_name}:"
        if rule_pattern not in content:
            logger.debug(f"Rule '{rule_name}' not found in Snakefile")
            return False

        # Extract rule metadata for validation
        metadata = get_rule_metadata(snakefile_path, rule_name)
        if not metadata:
            logger.warning(f"Rule '{rule_name}' found but metadata extraction failed")
            return False

        # Basic structure validation
        required_sections = ['input', 'output']
        for section in required_sections:
            if section not in metadata:
                logger.warning(f"Rule '{rule_name}' missing required section: {section}")
                return False

        logger.debug(f"Rule '{rule_name}' validation successful")
        return True

    except Exception as e:
        logger.error(f"Error validating rule '{rule_name}': {e}")
        return False


def get_rule_metadata(snakefile_path: str, rule_name: str) -> Optional[Dict[str, Any]]:
    """
    Extract rule metadata including inputs, outputs, and parameters.

    Args:
        snakefile_path: Path to the Snakefile
        rule_name: Name of the rule to analyze

    Returns:
        Dict containing rule metadata or None if extraction fails

    Example:
        metadata = get_rule_metadata("/path/to/Snakefile", "Heatmap")
        # Returns: {
        #     "inputs": ["expression.csv", "trajectory.txt", "input.txt"],
        #     "outputs": ["Heatmap.json"],
        #     "parameters": ["Top_Genes", "Sample_Size"],
        #     "shell_command": "Rscript ..."
        # }
    """
    try:
        with open(snakefile_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        rule_pattern = f"rule {rule_name}:"
        metadata = {
            'inputs': [],
            'outputs': [],
            'parameters': [],
            'shell_command': None
        }

        is_target_rule = False
        current_section = None
        current_indent = 0

        for line in lines:
            if line.startswith(rule_pattern):
                is_target_rule = True
                current_indent = len(line) - len(line.lstrip())
                continue

            if not is_target_rule:
                continue

            if line.strip() == "":
                continue

            line_indent = len(line) - len(line.lstrip())

            # End of rule block
            if line_indent <= current_indent and line.strip():
                break

            stripped_line = line.strip()

            # Identify sections
            if stripped_line.endswith(':'):
                section_name = stripped_line[:-1]
                if section_name in ['input', 'output', 'params', 'shell']:
                    current_section = section_name
                continue

            # Parse section content
            if current_section == 'input':
                # Extract input file patterns
                if '=' in stripped_line:
                    file_pattern = stripped_line.split('=')[1].strip().strip('"')
                    # Extract filename from pattern
                    filename = _extract_filename_from_pattern(file_pattern)
                    if filename:
                        metadata['inputs'].append(filename)

            elif current_section == 'output':
                # Extract output file patterns
                if '=' in stripped_line:
                    file_pattern = stripped_line.split('=')[1].strip().strip('"')
                    filename = _extract_filename_from_pattern(file_pattern)
                    if filename:
                        metadata['outputs'].append(filename)

            elif current_section == 'params':
                # Extract parameter names
                if '=' in stripped_line:
                    param_name = stripped_line.split('=')[0].strip()
                    metadata['parameters'].append(param_name)

            elif current_section == 'shell':
                # Extract shell command
                if stripped_line.startswith('"') and stripped_line.endswith('"'):
                    metadata['shell_command'] = stripped_line.strip('"')

        logger.debug(f"Extracted metadata for rule '{rule_name}': {metadata}")
        return metadata

    except Exception as e:
        logger.error(f"Error extracting metadata for rule '{rule_name}': {e}")
        return None


def _get_available_rules(snakefile_path: str) -> List[str]:
    """Extract list of available rule names from a Snakefile."""
    try:
        with open(snakefile_path, 'r', encoding='utf-8') as file:
            content = file.read()

        # Find all rule definitions
        rule_pattern = r'^rule\s+(\w+):'
        matches = re.findall(rule_pattern, content, re.MULTILINE)
        return matches

    except Exception as e:
        logger.error(f"Error getting available rules: {e}")
        return []


def _apply_replacements(content: str, parameters: Optional[Dict[str, Any]],
                       file_mappings: Optional[Dict[str, str]]) -> str:
    """
    Apply parameter and file mapping replacements to extracted rule content.

    Args:
        content: The Snakefile rule content to process
        parameters: Dict of parameter name/value pairs to replace
        file_mappings: Dict mapping file placeholders to absolute paths
                      e.g., {"input.txt": "/absolute/path/to/expression.csv"}
                      These replace placeholders like {input.txt} with actual file paths

    Returns:
        str: Content with all placeholders replaced

    Raises:
        SnakefileGenerationError: If placeholder replacement fails
    """
    try:
        # Apply parameter replacements
        if parameters:
            for param_name, param_value in parameters.items():
                # Replace parameter placeholders like {Top Genes} or {Sample Size}
                # Handle various naming conventions - covers all cases without duplication
                placeholder_patterns = [
                    f"{{{param_name}}}",  # Exact match
                    f"{{{param_name.replace('_', ' ')}}}",  # Underscore to space
                    f"{{{param_name.replace(' ', '_')}}}",  # Space to underscore
                ]

                for pattern in placeholder_patterns:
                    if pattern in content:
                        content = content.replace(pattern, str(param_value))
                        logger.debug(f"Replaced parameter placeholder: {pattern} -> {param_value}")

        # Apply file mapping replacements
        if file_mappings:
            for placeholder, actual_path in file_mappings.items():
                # Replace file placeholders in input/output sections
                placeholder_patterns = [
                    f"{{{placeholder}}}",  # Standard placeholder format
                    f"/{{{placeholder}}}",  # With leading slash
                    f'"{{{placeholder}}}"',  # In quotes
                    f"'{{{placeholder}}}'",  # In single quotes
                ]

                for pattern in placeholder_patterns:
                    if pattern in content:
                        # Preserve quotes if they exist
                        if pattern.startswith('"') and pattern.endswith('"'):
                            replacement = f'"{actual_path}"'
                        elif pattern.startswith("'") and pattern.endswith("'"):
                            replacement = f"'{actual_path}'"
                        else:
                            replacement = actual_path

                        content = content.replace(pattern, replacement)
                        logger.debug(f"Replaced file placeholder: {pattern} -> {replacement}")

        # Check for any remaining placeholders that weren't replaced
        remaining_placeholders = re.findall(r'\{[^}]+\}', content)
        if remaining_placeholders:
            logger.warning(f"Unreplaced placeholders found: {remaining_placeholders}")
            logger.debug(f"Available parameters: {list(parameters.keys()) if parameters else []}")
            logger.debug(f"Available file mappings: {list(file_mappings.keys()) if file_mappings else []}")

        return content

    except Exception as e:
        logger.error(f"Error applying replacements: {e}")
        # Re-raise the exception to prevent silent failures
        raise SnakefileGenerationError(
            plugin_name="unknown",
            script_name="unknown",
            error_details=f"Placeholder replacement failed: {str(e)}"
        ) from e


def _extract_filename_from_pattern(file_pattern: str) -> Optional[str]:
    """Extract filename from Snakemake file pattern."""
    try:
        # Handle patterns like "user/{user_name}/data/{filename}"
        if "/" in file_pattern:
            # Get the last part of the path
            filename = file_pattern.split("/")[-1]
        else:
            filename = file_pattern

        # Remove curly braces if present
        filename = filename.strip("{}")

        # Return None for template variables
        if filename in ['user_name', 'workflow_id', 'algorithm_id']:
            return None

        return filename

    except Exception:
        return None
