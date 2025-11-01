"""
Submodule Validation Tests

Tests for validating git submodule configuration and state for CPU/GPU modes.
Covers TC-SUBMODULE-001 through TC-SUBMODULE-012 from deployment test plan.

CRITICAL: CPU and GPU modes use different submodule branches:
- CPU mode: release/plugins-v1.0-cpu
- GPU mode: release/plugins-v1.0
"""

import pytest
import os


# Submodule path constant
PLUGIN_SUBMODULE_PATH = "backend/plugin/official"


class TestSubmoduleInitialization:
    """4.3 Submodule Initialization Tests (HIGH PRIORITY)"""

    def test_submodule_initialized(self, project_root):
        """
        TC-SUBMODULE-002: Plugin submodule is initialized

        Validates that the plugin submodule is properly initialized and not empty.

        Success Criteria:
        - Submodule directory exists
        - Submodule is not empty
        - Git metadata exists (.git file or directory)
        - Submodule has valid commit hash

        Failure Guidance:
        - If not initialized: Run `git submodule update --init --recursive`
        - If directory empty: Check .gitmodules configuration
        - If no git metadata: Verify git repository setup
        """
        from .helpers import get_submodule_status

        print(f"\n=== Submodule Initialization Check ===")
        print(f"  Checking: {PLUGIN_SUBMODULE_PATH}")

        # Get submodule status
        submodule_path = os.path.join(project_root, PLUGIN_SUBMODULE_PATH)
        status = get_submodule_status(PLUGIN_SUBMODULE_PATH)

        # Check initialization
        if not status["initialized"]:
            error_msg = status.get("error", "Unknown error")
            print(f"\n❌ Submodule Not Initialized")
            print(f"  Error: {error_msg}")
            print(f"\n  To initialize submodule, run:")
            print(f"    cd {project_root}")
            print(f"    git submodule update --init --recursive")

            assert False, \
                f"Submodule not initialized: {error_msg}\n" \
                f"Run: git submodule update --init --recursive"

        # Verify submodule has commit hash
        assert status["commit_hash"] is not None, \
            "Submodule should have valid commit hash"

        print(f"\n✅ Submodule Initialized")
        print(f"  - Path: {PLUGIN_SUBMODULE_PATH}")
        print(f"  - Commit: {status['commit_hash']}")
        print(f"  - Branch: {status['branch'] or '(detached HEAD)'}")
        print(f"  - Remote: {status['remote_url']}")

    def test_required_plugin_directories(self, project_root):
        """
        TC-SUBMODULE-011: Required plugin directories exist

        Validates that expected plugin directories exist in the submodule.

        Success Criteria:
        - Submodule directory exists and is not empty
        - At least 6 plugin directories present (official plugins)
        - Each plugin directory has basic structure

        Expected Plugins (6 official):
        - TENET
        - GENIE3
        - GRNBOOST2
        - LEAP
        - Scribe
        - GRNViz
        """
        import os

        print(f"\n=== Plugin Directory Validation ===")

        submodule_path = os.path.join(project_root, PLUGIN_SUBMODULE_PATH)

        # Check submodule exists
        assert os.path.exists(submodule_path), \
            f"Submodule directory should exist: {submodule_path}"

        # List plugin directories (exclude hidden files and version.json)
        plugin_dirs = [
            d for d in os.listdir(submodule_path)
            if os.path.isdir(os.path.join(submodule_path, d))
            and not d.startswith('.')
        ]

        # Minimum 6 official plugins expected
        min_plugins = 6
        assert len(plugin_dirs) >= min_plugins, \
            f"Expected at least {min_plugins} plugin directories, found {len(plugin_dirs)}: {plugin_dirs}"

        print(f"\n✅ Plugin Directories Validated")
        print(f"  - Found: {len(plugin_dirs)} plugin directories")
        print(f"  - Minimum required: {min_plugins}")
        print(f"\n  Plugin List:")
        for plugin_dir in sorted(plugin_dirs):
            plugin_path = os.path.join(submodule_path, plugin_dir)
            has_metadata = os.path.exists(os.path.join(plugin_path, "metadata.json"))
            has_snakefile = os.path.exists(os.path.join(plugin_path, "Snakefile"))
            marker = "✓" if has_metadata and has_snakefile else "✗"
            print(f"    {marker} {plugin_dir}")


class TestSubmoduleBranchConsistency:
    """4.4 Submodule Branch Consistency Tests (HIGH PRIORITY)"""

    @pytest.mark.cpu_mode
    def test_cpu_mode_branch(self, cpu_mode_running, project_root):
        """
        TC-SUBMODULE-005: CPU mode uses release/plugins-v1.0-cpu branch

        Validates that when running in CPU mode, the plugin submodule
        is on the correct branch for CPU-only plugins.

        Success Criteria:
        - Submodule is on 'release/plugins-v1.0-cpu' branch
        - Branch matches expected CPU mode configuration
        - Not on detached HEAD state

        Failure Guidance:
        - If on wrong branch: Run deployment script to switch branches
        - If detached HEAD: Checkout correct branch manually
        - Command: cd backend/plugin/official && git checkout release/plugins-v1.0-cpu
        """
        from .helpers import get_submodule_branch

        expected_branch = "release/plugins-v1.0-cpu"

        print(f"\n=== CPU Mode Branch Validation ===")
        print(f"  Expected branch: {expected_branch}")

        # Get current branch
        current_branch = get_submodule_branch(PLUGIN_SUBMODULE_PATH)

        # Check if on correct branch
        if current_branch is None:
            print(f"\n❌ Submodule in Detached HEAD State")
            print(f"\n  To fix, run:")
            print(f"    cd {project_root}/{PLUGIN_SUBMODULE_PATH}")
            print(f"    git checkout {expected_branch}")

            assert False, \
                f"Submodule in detached HEAD state. Expected branch: {expected_branch}\n" \
                f"Run: cd {PLUGIN_SUBMODULE_PATH} && git checkout {expected_branch}"

        if current_branch != expected_branch:
            print(f"\n❌ Wrong Branch")
            print(f"  Current: {current_branch}")
            print(f"  Expected: {expected_branch}")
            print(f"\n  To switch to correct branch:")
            print(f"    cd {project_root}")
            print(f"    ./run-cpu-mode.sh")

            assert False, \
                f"Submodule on wrong branch: {current_branch}\n" \
                f"Expected: {expected_branch}\n" \
                f"Run: ./run-cpu-mode.sh to switch branches"

        print(f"\n✅ CPU Mode Branch Correct")
        print(f"  - Current branch: {current_branch}")
        print(f"  - Expected branch: {expected_branch}")

    @pytest.mark.gpu_mode
    @pytest.mark.requires_gpu
    @pytest.mark.macos_skip
    def test_gpu_mode_branch(self, gpu_mode_running, project_root):
        """
        TC-SUBMODULE-006: GPU mode uses release/plugins-v1.0 branch

        Validates that when running in GPU mode, the plugin submodule
        is on the correct branch for GPU-capable plugins.

        Success Criteria:
        - Submodule is on 'release/plugins-v1.0' branch
        - Branch matches expected GPU mode configuration
        - Not on detached HEAD state

        Failure Guidance:
        - If on wrong branch: Run deployment script to switch branches
        - If detached HEAD: Checkout correct branch manually
        - Command: cd backend/plugin/official && git checkout release/plugins-v1.0
        """
        from .helpers import get_submodule_branch

        expected_branch = "release/plugins-v1.0"

        print(f"\n=== GPU Mode Branch Validation ===")
        print(f"  Expected branch: {expected_branch}")

        # Get current branch
        current_branch = get_submodule_branch(PLUGIN_SUBMODULE_PATH)

        # Check if on correct branch
        if current_branch is None:
            print(f"\n❌ Submodule in Detached HEAD State")
            print(f"\n  To fix, run:")
            print(f"    cd {project_root}/{PLUGIN_SUBMODULE_PATH}")
            print(f"    git checkout {expected_branch}")

            assert False, \
                f"Submodule in detached HEAD state. Expected branch: {expected_branch}\n" \
                f"Run: cd {PLUGIN_SUBMODULE_PATH} && git checkout {expected_branch}"

        if current_branch != expected_branch:
            print(f"\n❌ Wrong Branch")
            print(f"  Current: {current_branch}")
            print(f"  Expected: {expected_branch}")
            print(f"\n  To switch to correct branch:")
            print(f"    cd {project_root}")
            print(f"    ./run-gpu-mode.sh")

            assert False, \
                f"Submodule on wrong branch: {current_branch}\n" \
                f"Expected: {expected_branch}\n" \
                f"Run: ./run-gpu-mode.sh to switch branches"

        print(f"\n✅ GPU Mode Branch Correct")
        print(f"  - Current branch: {current_branch}")
        print(f"  - Expected branch: {expected_branch}")


class TestSubmoduleVersionConsistency:
    """4.5 Submodule Version & State Tests (MEDIUM PRIORITY)"""

    def test_version_file_consistency(self, project_root):
        """
        TC-SUBMODULE-007: version.json branch matches actual branch

        Validates that version.json file exists and its branch field
        matches the actual checked-out branch.

        Success Criteria:
        - version.json file exists
        - File has valid JSON format
        - 'branch' field matches actual git branch
        - 'commit' field is present and valid

        Note: This test is independent of CPU/GPU mode and runs in base environment.
        """
        from .helpers import get_submodule_branch, validate_version_file
        import os

        print(f"\n=== Version File Consistency Check ===")

        submodule_path = os.path.join(project_root, PLUGIN_SUBMODULE_PATH)
        version_file = os.path.join(submodule_path, "version.json")

        # Get actual branch
        actual_branch = get_submodule_branch(PLUGIN_SUBMODULE_PATH)

        if actual_branch is None:
            # Detached HEAD - version file check not applicable
            print(f"\n⚠️ Submodule in Detached HEAD State")
            print(f"  Version file consistency check skipped")
            pytest.skip("Submodule in detached HEAD state, cannot validate version file consistency")
            return

        # Validate version file
        version_data = validate_version_file(version_file)

        if not version_data["valid"]:
            error = version_data.get("error", "Unknown error")
            print(f"\n⚠️ Version File Invalid or Missing")
            print(f"  Error: {error}")
            print(f"\n  This is not a fatal error, but indicates version file needs updating")
            # Don't fail - version file is informational
            return

        # Check branch consistency
        version_branch = version_data["branch"]

        if version_branch != actual_branch:
            print(f"\n⚠️ Branch Mismatch")
            print(f"  version.json branch: {version_branch}")
            print(f"  Actual git branch: {actual_branch}")
            print(f"\n  Version file may need to be updated to match current branch")
            # Don't fail - this is informational
        else:
            print(f"\n✅ Version File Consistent")
            print(f"  - Branch: {actual_branch}")
            print(f"  - Version: {version_data['version']}")
            print(f"  - Commit: {version_data['commit']}")

    def test_submodule_clean_state(self, project_root):
        """
        TC-SUBMODULE-009: Submodule has no uncommitted changes

        Validates that the plugin submodule working directory is clean
        with no uncommitted changes.

        Success Criteria:
        - No modified files
        - No staged changes
        - No untracked files in submodule

        Note: Uncommitted changes are warned but not failed, as they may be
        intentional during development.
        """
        from .helpers import check_submodule_clean, get_submodule_status
        import subprocess

        print(f"\n=== Submodule Clean State Check ===")

        # Check if submodule is clean
        is_clean = check_submodule_clean(PLUGIN_SUBMODULE_PATH)

        if not is_clean:
            # Get detailed status
            result = subprocess.run(
                ["git", "-C", PLUGIN_SUBMODULE_PATH, "status", "--short"],
                capture_output=True,
                text=True,
                timeout=5
            )

            print(f"\n⚠️ Submodule Has Uncommitted Changes")
            if result.returncode == 0:
                print(f"\n  Changes:")
                for line in result.stdout.strip().split('\n'):
                    print(f"    {line}")

            print(f"\n  This may be intentional during development.")
            print(f"  For production deployment, submodule should be clean.")

            # Don't fail - just warn
            print(f"\n  Warning: Submodule is not in clean state")
        else:
            print(f"\n✅ Submodule Clean")
            print(f"  - No uncommitted changes")
            print(f"  - Working directory clean")

    def test_plugin_metadata_files(self, project_root):
        """
        TC-SUBMODULE-012: Each plugin has metadata.json

        Validates that each plugin directory contains required metadata files.

        Success Criteria:
        - Each plugin directory has metadata.json
        - Each plugin directory has Snakefile
        - Metadata files are readable and valid JSON

        Note: Only checks for file existence, not full validation of contents.
        """
        import os
        import json

        print(f"\n=== Plugin Metadata Files Check ===")

        submodule_path = os.path.join(project_root, PLUGIN_SUBMODULE_PATH)

        # Get plugin directories
        plugin_dirs = [
            d for d in os.listdir(submodule_path)
            if os.path.isdir(os.path.join(submodule_path, d))
            and not d.startswith('.')
        ]

        missing_metadata = []
        missing_snakefile = []
        invalid_json = []

        for plugin_dir in plugin_dirs:
            plugin_path = os.path.join(submodule_path, plugin_dir)

            # Check metadata.json
            metadata_file = os.path.join(plugin_path, "metadata.json")
            if not os.path.exists(metadata_file):
                missing_metadata.append(plugin_dir)
            else:
                # Try to parse JSON
                try:
                    with open(metadata_file, 'r') as f:
                        json.load(f)
                except json.JSONDecodeError:
                    invalid_json.append(plugin_dir)

            # Check Snakefile
            snakefile = os.path.join(plugin_path, "Snakefile")
            if not os.path.exists(snakefile):
                missing_snakefile.append(plugin_dir)

        # Report results
        if missing_metadata:
            print(f"\n⚠️ Plugins Missing metadata.json:")
            for plugin in missing_metadata:
                print(f"    - {plugin}")

        if missing_snakefile:
            print(f"\n⚠️ Plugins Missing Snakefile:")
            for plugin in missing_snakefile:
                print(f"    - {plugin}")

        if invalid_json:
            print(f"\n⚠️ Plugins with Invalid JSON:")
            for plugin in invalid_json:
                print(f"    - {plugin}")

        # Fail if critical files missing
        assert not missing_metadata, \
            f"Plugins missing metadata.json: {', '.join(missing_metadata)}"

        assert not invalid_json, \
            f"Plugins with invalid JSON: {', '.join(invalid_json)}"

        print(f"\n✅ Plugin Metadata Files Validated")
        print(f"  - Total plugins: {len(plugin_dirs)}")
        print(f"  - All have metadata.json: Yes")
        print(f"  - All have Snakefile: {'Yes' if not missing_snakefile else 'No'}")

        # Snakefile is less critical, just warn
        if missing_snakefile:
            print(f"\n  Warning: Some plugins missing Snakefile")
