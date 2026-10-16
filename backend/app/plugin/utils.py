"""Backward-compat facade for the plugin utility helpers (PR-9 split).

The implementations moved into focused modules:

- ``plugin/files.py``    — path resolution & folder/file scaffolding
- ``plugin/metadata.py`` — dependency verification & drawflow templates
- ``plugin/runtime.py``  — Snakefile & Dockerfile generation
- ``plugin/builder.py``  — Docker image build/tag/check operations

Existing call sites and test patch targets that reference these names on
``app.plugin.utils`` keep working via the re-exports below. Docker build/tag
symbols now live on ``app.plugin.builder`` and are imported from there directly
(kept out of this facade to avoid a plugin<->worker import cycle).
"""
from app.plugin.files import (  # noqa: F401
    OFFICIAL_PLUGINS_DIR,
    LOCAL_PLUGINS_DIR,
    sanitize_plugin_name,
    get_plugin_path,
    list_available_plugins,
    resolve_plugin_file_path,
    is_plugin_editable,
    ensure_local_plugins_dir,
    create_plugin_folder,
    create_dependency_folder,
    create_reference_folder,
    get_reference_folders_list,
    get_reference_folder,
    create_metadata_file,
)
from app.plugin.metadata import (  # noqa: F401
    normalize_pkg_name,
    check_requirements_txt,
    check_environment_yml,
    check_renv_lock,
    verify_dependencies,
    generate_merged_plugin_drawflow,
    generate_plugin_drawflow_template,
)
from app.plugin.runtime import (  # noqa: F401
    normalize_param_name,
    normalize_string,
    generate_snakemake_code,
    extract_python_version_from_environment_yml,
    extract_r_version_from_renv_lock,
    generate_base_image_section,
    generate_env_setup_section,
    generate_system_packages_section,
    generate_dependency_copy_section,
    generate_micromamba_install_section,
    generate_env_variables_section,
    generate_micromamba_setup_section,
    generate_python_env_section,
    generate_snakemake_install_section,
    generate_python_packages_section,
    generate_r_installation_section,
    generate_workspace_setup_section,
    generate_copy_files_section,
    generate_entrypoint_section,
    generate_healthcheck_section,
    generate_cmd_section,
    generate_plugin_dockerfile,
)
