"""Plugin runtime/build artifact generation (split from ``plugin/utils.py`` in PR-9).

Generates the files needed to build & run a plugin:

- ``generate_snakemake_code``: the Snakefile that drives plugin execution
  (plus its ``normalize_param_name`` / ``normalize_string`` helpers).
- ``generate_plugin_dockerfile`` and the ``generate_*_section`` builders that
  assemble the plugin's Docker build recipe, with ``extract_*_version`` reading
  the pinned Python/R versions out of the dependency manifests.

Function signatures and behavior are unchanged from the original ``utils.py``.
"""
import os
import re
import json

import yaml
from fastapi import HTTPException


def normalize_param_name(name: str) -> str:
    """파라미터 이름에서 특수문자를 제거하고 언더스코어로 변환합니다."""
    # 공백을 언더스코어로 변환
    name = name.replace(' ', '_')
    # 알파벳, 숫자, 언더스코어만 허용하고 나머지는 제거
    name = re.sub(r'[^a-zA-Z0-9_]', '', name)
    return name

def normalize_string(text: str) -> str:
    """문자열에서 특수문자를 제거합니다."""
    # 알파벳, 숫자, 공백, 점(.), 언더스코어(_), 하이픈(-) 만 허용
    return re.sub(r'[^a-zA-Z0-9\s._-]', '', text)

def generate_snakemake_code(rules_data, output_folder_path, plugin_type=None, workflow_id=None, algo_id=None, viz_id=None):
    """
    Generate unified Snakemake code for both analysis and visualization plugins.
    
    Parameters:
        rules_data (dict): Rules data
        output_folder_path (str): Output folder path
        plugin_type (str): Plugin type ('ANALYSIS' or 'VISUALIZATION')
        workflow_id (str): Workflow ID for visualization plugins
        algo_id (str): Algorithm ID for input file resolution
        viz_id (str): Visualization node ID for output path
    
    Raises:
        HTTPException: If there's an error generating the Snakemake code
    """
    try:
        if not rules_data:
            raise ValueError("No rules data provided")

        snakemake_path = os.path.join(output_folder_path, "Snakefile")
        snakemake_code = "import os\n\n"  # Add import statement at the top
        
        # Configure paths based on plugin type
        # When workflow_id, algo_id, viz_id are provided, use them; otherwise use placeholders
        if plugin_type == 'VISUALIZATION':
            # For visualization plugins: input from algorithm results, output to visualization directory
            if workflow_id and algo_id and viz_id:
                # Specific values provided - for runtime execution
                input_path = f"user/{{user_name}}/workflow_{workflow_id}/algorithm_{algo_id}/results"
                output_path = f"user/{{user_name}}/workflow_{workflow_id}/visualization_{viz_id}/results"
                logs_path = f"user/{{user_name}}/workflow_{workflow_id}/visualization_{viz_id}/logs"
            else:
                # No specific values - generate template with placeholders
                input_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/results"
                output_path = "user/{user_name}/workflow_{workflow_id}/visualization_{visualization_id}/results"
                logs_path = "user/{user_name}/workflow_{workflow_id}/visualization_{visualization_id}/logs"
            unique_input_path = input_path  # Visualization plugins don't use unique inputs from data folder
        else:
            # For analysis plugins: input from data, output to algorithm results
            if workflow_id and algo_id:
                # Specific values provided - for runtime execution
                input_path = f"user/{{user_name}}/workflow_{workflow_id}/algorithm_{algo_id}/results"
                output_path = f"user/{{user_name}}/workflow_{workflow_id}/algorithm_{algo_id}/results"
                logs_path = f"user/{{user_name}}/workflow_{workflow_id}/algorithm_{algo_id}/logs"
            else:
                # No specific values - generate template with placeholders
                input_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/results"
                output_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/results"
                logs_path = "user/{user_name}/workflow_{workflow_id}/algorithm_{algorithm_id}/logs"
            unique_input_path = "user/{user_name}/data"

        # Find all outputs across all rules to determine unique inputs not present in outputs
        all_outputs = {out for rule in rules_data.values() for out in rule['output']}
        all_inputs = {inp for rule in rules_data.values() for inp in rule['input']}
        unique_inputs = [inp for inp in all_inputs if inp not in all_outputs]
        scripts_log_word = "> {log.stdout} 2> {log.stderr}"

        # Iterate through each rule in the dictionary
        for rule_id, rule in rules_data.items():
            # All rules go into the same unified Snakefile
            print(f"Adding rule to Snakefile: {rule['name']} (plugin_type: {plugin_type})")
            
            # Start rule block with rule name
            snakemake_code += f"rule {normalize_param_name(rule['name'])}:\n"

            # Input section with optional input handling
            if 'input' in rule and rule['input']:
                def get_input_param_name(inp, rule_params):
                    for param in rule_params:
                        if param['type'] == 'inputFile' and param.get('defaultValue') == inp:
                            return normalize_param_name(param['name'])
                    return normalize_param_name(os.path.splitext(os.path.basename(inp))[0])

                input_files = ",\n        ".join([
                    f"{get_input_param_name(inp, rule.get('parameters', []))}=\"{input_path}/{{{inp}}}\"" 
                    if 'target' in get_input_param_name(inp, rule.get('parameters', [])) else
                    f"{get_input_param_name(inp, rule.get('parameters', []))}=\"{unique_input_path}/{{{inp}}}\"" 
                    if inp in unique_inputs else
                    f"{get_input_param_name(inp, rule.get('parameters', []))}=\"{input_path}/{inp}\""
                    for inp in rule['input']
                    if not any(param['type'] == 'optionalInputFile' and param.get('defaultValue') == inp 
                              for param in rule.get('parameters', []))
                ])
                if input_files:  # Only add input section if there are non-optional input files
                    snakemake_code += f"    input:\n        {input_files}\n"

            # Output section
            if 'output' in rule and rule['output']:
                def get_output_param_name(out, rule_params):
                    for param in rule_params:
                        if param['type'] == 'outputFile' and param.get('defaultValue') == out:
                            return normalize_param_name(param['name'])
                    return normalize_param_name(os.path.splitext(os.path.basename(out))[0])

                # Use appropriate output path based on plugin type
                output_files = ",\n        ".join([
                    f"{get_output_param_name(out, rule.get('parameters', []))}=\"{output_path}/{out}\""
                    for out in rule['output']
                ])
                snakemake_code += f"    output:\n        {output_files}\n"

            # Params section
            if 'parameters' in rule and rule['parameters']:
                param_list = []
                for param in rule['parameters']:
                    normalized_name = normalize_param_name(param['name'])
                    if param['name'] == "clusters" and param['type'] == "h5adParameter":
                        param_list.append(f'clusters=lambda wc: ";".join({{{param["name"]}}})')
                    elif param['name'] == "ScatterPlot" and param['type'] == 'string':
                        param_list.append(
                            'ScatterPlot=lambda wildcards: {ScatterPlot} if os.path.exists({ScatterPlot}) else "None"'
                        )
                    elif "UMAP" in param['name'] and param['type'] == 'h5adParameter':
                        param_list.append(
                            'UMAP_lasso=lambda wildcards: {UMAP lasso} if os.path.exists({UMAP lasso}) else "None"'
                        )
                    elif param['type'] == 'optionalInputFile':
                        param_list.append(
                            f"{normalized_name}=lambda wildcards: \"{unique_input_path}/{{{param['defaultValue']}}}\" "
                            f"if os.path.exists(\"{unique_input_path}/{{{param['defaultValue']}}}\") else \"None\""
                        )
                    elif param['type'] != 'inputFile' and param['type'] != 'outputFile':
                        param_list.append(f"{normalized_name}={{{param['name']}}}")

                if param_list:
                    param_list_str = ",\n        ".join(param_list)
                    snakemake_code += f"    params:\n        {param_list_str}\n"

            # Log section
            snakemake_code += f"    log:\n"
            snakemake_code += f"        stdout=\"{logs_path}/{normalize_param_name(rule['name'])}.stdout\",\n"
            snakemake_code += f"        stderr=\"{logs_path}/{normalize_param_name(rule['name'])}.stderr\"\n"

            # Shell section based on script type
            if 'script' in rule and rule['script']:
                # Script path 설정 (Docker 내부 경로)
                script_path = f"/scripts/{normalize_string(rule['script'])}"

                # Shell command 설정 (Python/R 분기)
                if rule['script'].endswith('.py'):
                    shell_command = f"/opt/micromamba/envs/plugin_env/bin/python {script_path}"
                elif rule['script'].endswith('.R'):
                    # Simplified R script execution using micromamba environment
                    shell_command = f"export R_HOME=/opt/micromamba/envs/plugin_env/lib/R && " \
                                  f"export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE && " \
                                  f"/opt/micromamba/envs/plugin_env/bin/Rscript {script_path}"
                else:
                    shell_command = script_path

                # Add parameters to shell command with input/output/params distinction
                if 'parameters' in rule and rule['parameters']:
                    param_list = []
                    for param in rule['parameters']:
                        normalized_name = normalize_param_name(param['name'])
                        if param['type'] == 'inputFile':
                            param_list.append(f"{{input.{normalized_name}}}")
                        elif param['type'] == 'optionalInputFile':
                            param_list.append(f"{{params.{normalized_name}}}")
                        elif param['type'] == 'outputFile':
                            param_list.append(f"{{output.{normalized_name}}}")
                        elif param['name'] == "clusters" and param['type'] == "h5adParameter":
                            param_list.append(f"'{{{normalize_string(f'params.{normalized_name}')}}}'")
                        else:
                            param_list.append(f"{{{normalize_string(f'params.{normalized_name}')}}}")

                    param_list_str = " ".join(param_list)
                    shell_command = f"{shell_command} {param_list_str}"

                snakemake_code += f"    shell:\n        \"{shell_command} {scripts_log_word}\"\n"

            snakemake_code += "\n"  # Add a newline for separation between rules

        # Write the unified Snakefile
        with open(snakemake_path, 'w') as file:
            file.write(snakemake_code)
            
        print(f"Generated unified Snakefile at: {snakemake_path} (plugin_type: {plugin_type})")

    except ValueError as ve:
        raise HTTPException(status_code=400, detail=str(ve))
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to generate Snakemake code: {str(e)}"
        )


def extract_python_version_from_environment_yml(env_path: str) -> str:
    """environment.yml 파일에서 Python 버전을 추출합니다."""
    if not os.path.isfile(env_path):
        return None
    with open(env_path, 'r') as f:
        data = yaml.safe_load(f)
    dependencies = data.get('dependencies', [])
    for dep in dependencies:
        if isinstance(dep, str) and dep.startswith("python="):
            return dep.split("=")[1]
    return None

def extract_r_version_from_renv_lock(renv_path: str) -> str:
    """renv.lock 파일에서 R 버전을 추출합니다."""
    if not os.path.isfile(renv_path):
        return None
    with open(renv_path, 'r') as f:
        data = json.load(f)
    r_info = data.get("R", {})
    if isinstance(r_info, dict):
        return r_info.get("Version", None)
    return None


def generate_base_image_section(use_gpu=True):
    """Base image section generation"""
    lines = []
    if use_gpu:
        lines.append("FROM --platform=linux/amd64 nvidia/cuda:12.1.0-cudnn8-devel-ubuntu20.04")
    else:
        lines.append("ARG TARGETPLATFORM=linux/amd64")
        lines.append("ARG BUILDPLATFORM=linux/amd64")
        lines.append("")
        lines.append("FROM --platform=$TARGETPLATFORM debian:bullseye-slim")
    return lines


def generate_env_setup_section():
    """Environment variables setup section"""
    return [
        "",
        "ENV DEBIAN_FRONTEND=noninteractive",
        "ENV TZ=Asia/Seoul"
    ]


def generate_system_packages_section(use_gpu=True, has_r=False):
    """Enhanced system packages installation section based on Scribe requirements"""
    lines = [
        "",
        "RUN apt-get update && apt-get install -y \\"
    ]
    
    # Common packages
    common_packages = [
        "    build-essential gcc g++ gfortran make \\",
        "    libssl-dev libcurl4-openssl-dev libxml2-dev libxslt1-dev \\",
        "    libjpeg-dev libpng-dev libfreetype6-dev libtiff5-dev \\",
        "    libx11-dev libxt-dev" + (" xorg-dev \\" if has_r else " \\"),
        "    libglu1-mesa-dev \\",
        "    libharfbuzz-dev libfribidi-dev \\",
        "    libglpk-dev \\",
        "    curl wget unzip git \\",
        "    python3 python3-pip python3-venv \\",
        "    software-properties-common \\",
        "    gnupg ca-certificates \\",
        "    zlib1g-dev" + (" liblzma-dev \\" if has_r else "")
    ]
    
    lines.extend(common_packages)
    
    # Comprehensive R ecosystem packages (optimized for Scribe)
    if has_r:
        r_packages = [
            "    rsync \\",
            "    libv8-dev libnode-dev \\",
            "    libudunits2-dev libgdal-dev \\",
            "    libcairo2-dev \\",
            "    libfontconfig1-dev \\", 
            "    libgeos-dev \\",
            "    libproj-dev \\",
            "    libmagick++-dev \\",
            "    libpoppler-cpp-dev \\",
            "    librsvg2-dev \\"
        ]
        lines.extend(r_packages)
    
    # 마지막 패키지 라인이 백슬래시로 끝나지 않으면 백슬래시 추가
    if not lines[-1].endswith(" \\"):
        lines[-1] = lines[-1] + " \\"
    
    lines.append("    && apt-get clean && rm -rf /var/lib/apt/lists/*")
    return lines


def generate_dependency_copy_section():
    """Dependency folder copy section"""
    return [
        "",
        "COPY dependency/ /workspace/dependency/"
    ]


def generate_micromamba_install_section(use_gpu=False):
    """Micromamba installation section - dynamic platform-specific installation"""
    if use_gpu:
        # GPU 버전은 linux/amd64 고정
        return [
            "",
            "RUN mkdir -p /usr/local/bin && \\",
            "    cd /tmp && \\",
            "    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj && \\",
            "    cp bin/micromamba /usr/local/bin/micromamba && \\",
            "    chmod +x /usr/local/bin/micromamba && \\",
            "    rm -rf /tmp/bin /tmp/info && \\",
            "    /usr/local/bin/micromamba --version"
        ]
    else:
        # CPU 버전은 플랫폼별 동적 설치 (기본값 포함)
        return [
            "",
            "RUN PLATFORM=${TARGETPLATFORM:-linux/amd64} && \\",
            "    case $PLATFORM in \\",
            "        \"linux/amd64\")  MICROMAMBA_ARCH=linux-64  ;; \\",
            "        \"linux/arm64\")  MICROMAMBA_ARCH=linux-aarch64  ;; \\",
            "        \"linux/arm/v7\") MICROMAMBA_ARCH=linux-armv7l  ;; \\",
            "        *) echo \"Unsupported platform: $PLATFORM, defaulting to linux-64\" && MICROMAMBA_ARCH=linux-64 ;; \\",
            "    esac && \\",
            "    mkdir -p /usr/local/bin && \\",
            "    curl -Ls \"https://micro.mamba.pm/api/micromamba/${MICROMAMBA_ARCH}/latest\" | tar -xvj -C /tmp && \\",
            "    cp /tmp/bin/micromamba /usr/local/bin/micromamba && \\",
            "    chmod +x /usr/local/bin/micromamba && \\",
            "    rm -rf /tmp/bin /tmp/info && \\",
            "    /usr/local/bin/micromamba --version"
        ]


def generate_env_variables_section(has_r=False):
    """Environment variables configuration section"""
    lines = [
        "",
        "ENV MAMBA_ROOT_PREFIX=/opt/micromamba",
        "ENV MAMBA_EXE=/usr/local/bin/micromamba",
        "ENV PATH=/usr/local/bin:$MAMBA_ROOT_PREFIX/envs/plugin_env/bin:$MAMBA_ROOT_PREFIX/bin:$PATH"
    ]
    
    if has_r:
        lines.append("ENV R_HOME=/opt/micromamba/envs/plugin_env/lib/R")
    
    return lines


def generate_micromamba_setup_section():
    """Micromamba environment directory creation section"""
    return [
        "",
        "RUN mkdir -p $MAMBA_ROOT_PREFIX && \\",
        "    mkdir -p $MAMBA_ROOT_PREFIX/envs && \\",
        "    mkdir -p $MAMBA_ROOT_PREFIX/pkgs && \\",
        "    mkdir -p $MAMBA_ROOT_PREFIX/etc/profile.d"
    ]


def generate_python_env_section(python_version="3.10"):
    """Python environment creation section"""
    return [
        "",
        f"RUN /usr/local/bin/micromamba create -y -n plugin_env python={python_version} -c conda-forge --root-prefix $MAMBA_ROOT_PREFIX"
    ]


def generate_snakemake_install_section():
    """Snakemake and essential packages installation section"""
    return [
        "",
        "RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \\",
        "    pip install --no-cache-dir \\",
        "    'snakemake==7.14.0' \\",
        "    'pulp==2.7.0' \\",
        "    'tabulate==0.8.10'"
    ]


def generate_python_packages_section(has_requirements=False):
    """Python packages installation section"""
    if not has_requirements:
        return []
    
    return [
        "",
        "RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \\",
        "    pip install --no-cache-dir -r /workspace/dependency/requirements.txt || true"
    ]


def generate_r_installation_section(r_version="4.4.2", has_renv=False, use_gpu=False):
    """Enhanced R installation section - efficient package management with rsync"""
    lines = []
    
    # Comprehensive Micromamba R packages installation
    lines.extend([
        "",
        f"RUN /usr/local/bin/micromamba install -y -n plugin_env \\",
        f"    r-base={r_version} \\",
        "    r-essentials \\",
        "    r-renv \\",
        "    r-jsonlite \\",
        "    zlib \\",
        "    xz \\",
        "    libpng \\",
        "    libjpeg-turbo \\",
        "    libcurl \\",
        "    libxml2 \\",
        "    libxslt \\",
        "    cairo \\",
        "    pango \\",
        "    harfbuzz \\",
        "    fribidi \\",
        "    freetype \\",
        "    fontconfig \\",
        "    libtiff \\",
        "    pkg-config \\",
        "    -c conda-forge \\",
        "    -r $MAMBA_ROOT_PREFIX && \\",
        "    ln -sf $MAMBA_ROOT_PREFIX/envs/plugin_env/bin/R /usr/local/bin/R && \\",
        "    ln -sf $MAMBA_ROOT_PREFIX/envs/plugin_env/bin/Rscript /usr/local/bin/Rscript"
    ])
    
    # R configuration
    lines.extend([
        "",
        "RUN echo 'options(repos = c(CRAN = \"https://cloud.r-project.org\"), download.file.method = \"libcurl\")' > /root/.Rprofile && \\",
        "    /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \\",
        "    Rscript -e \"install.packages(c('renv', 'BiocManager'), dependencies = TRUE)\""
    ])
    
    # Enhanced renv packages installation
    if has_renv:
        lines.extend([
            "",
            "RUN if [ -f \"/workspace/dependency/renv.lock\" ]; then \\",
            "    echo 'Installing R packages from renv.lock...' && \\",
            "    /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \\",
            "    Rscript -e \"\\",
            "        # Set micromamba environment as primary library path \\",
            "        target_lib <- '/opt/micromamba/envs/plugin_env/lib/R/library'; \\",
            "        .libPaths(c(target_lib, .libPaths())); \\",
            "        cat('Library paths configured:', paste(.libPaths(), collapse=', '), '\\\\n'); \\",
            "        \\",
            "        # Change to dependency directory for renv operations \\",
            "        setwd('/workspace/dependency'); \\",
            "        \\",
            "        # Load renv and configure \\",
            "        library(renv); \\",
            "        options(renv.config.cache.enabled = FALSE); \\",
            "        options(renv.config.auto.snapshot = FALSE); \\",
            "        \\",
            "        # Initialize renv project \\",
            "        cat('Initializing renv project...\\\\n'); \\",
            "        renv::init(bare = TRUE, restart = FALSE); \\",
            "        \\",
            "        # Restore packages to renv library \\",
            "        cat('Restoring packages from renv.lock...\\\\n'); \\",
            "        tryCatch({ \\",
            "            renv::restore(lockfile = '/workspace/dependency/renv.lock', confirm = FALSE); \\",
            "            cat('✅ R packages successfully installed to renv library\\\\n'); \\",
            "        }, error = function(e) { \\",
            "            cat('⚠️ Some packages may have failed, but continuing...\\\\n'); \\",
            "            cat('Error details:', conditionMessage(e), '\\\\n'); \\",
            "        }); \\",
            "    \" && \\",
            "    echo 'renv installation phase completed' ; \\",
            "fi",
            "",
            "RUN if [ -d \"/workspace/dependency/renv/library\" ]; then \\",
            "    echo '=== Stage 1: Starting complete rsync copy ===' && \\",
            "    RENV_LIB_PATH=\"/workspace/dependency/renv/library/linux-debian-bullseye/R-4.4/x86_64-conda-linux-gnu\" && \\",
            "    TARGET_LIB_PATH=\"/opt/micromamba/envs/plugin_env/lib/R/library\" && \\",
            "    if [ -d \"$RENV_LIB_PATH\" ]; then \\",
            "        echo \"Copy source: $RENV_LIB_PATH\" && \\",
            "        echo \"Copy target: $TARGET_LIB_PATH\" && \\",
            "        echo \"Packages to copy:\" && \\",
            "        ls -1 \"$RENV_LIB_PATH\" | wc -l && echo \" packages\" && \\",
            "        echo \"Starting rsync copy...\" && \\",
            "        rsync -av --progress \"$RENV_LIB_PATH/\" \"$TARGET_LIB_PATH/\" && \\",
            "        echo \"✅ rsync copy completed\" && \\",
            "        echo \"Post-copy package verification:\" && \\",
            "        ls -1 \"$TARGET_LIB_PATH\" | wc -l && echo \" packages\" ; \\",
            "    else \\",
            "        echo \"⚠️ Cannot find renv library path: $RENV_LIB_PATH\" ; \\",
            "    fi ; \\",
            "else \\",
            "    echo \"⚠️ renv library directory does not exist\" ;"
            "fi",
            "",
            "RUN echo '=== Stage 2: Cleanup temporary files for space optimization ===' && \\",
            "    echo '=== Stage 2: Starting temporary file cleanup ===' && \\",
            "    # Remove renv temporary library directory \\",
            "    if [ -d \"/workspace/dependency/renv/library\" ]; then \\",
            "        echo 'Removing renv temporary library...' && \\",
            "        RENV_LIB_SIZE=$(du -sh /workspace/dependency/renv/library 2>/dev/null | cut -f1 || echo \"unknown\") && \\",
            "        echo \"renv library size to remove: $RENV_LIB_SIZE\" && \\",
            "        rm -rf /workspace/dependency/renv/library && \\",
            "        echo '✅ renv temporary library removal completed' ; \\",
            "    fi && \\",
            "    # Remove renv cache \\",
            "    if [ -d \"/root/.cache/R/renv\" ]; then \\",
            "        echo 'Removing renv cache...' && \\",
            "        CACHE_SIZE=$(du -sh /root/.cache/R/renv 2>/dev/null | cut -f1 || echo \"unknown\") && \\",
            "        echo \"Cache size to remove: $CACHE_SIZE\" && \\",
            "        rm -rf /root/.cache/R/renv && \\",
            "        echo '✅ renv cache removal completed' ; \\",
            "    fi && \\",
            "    # Remove rsync package \\",
            "    echo 'Removing rsync package...' && \\",
            "    apt-get remove -y rsync && \\",
            "    apt-get autoremove -y && \\",
            "    echo '✅ rsync package removal completed' && \\",
            "    echo '=== Stage 2: Temporary file cleanup completed ==='"
        ])
    
    lines.extend([
        "",
        "RUN /usr/local/bin/micromamba run -n plugin_env -r $MAMBA_ROOT_PREFIX \\",
        "    Rscript -e \"if (!requireNamespace('optparse', quietly = TRUE)) install.packages('optparse', dependencies = TRUE)\"",
        "",
        "ENV RENV_CONFIG_AUTOLOADER_ENABLED=FALSE"
    ])
    
    return lines

def generate_workspace_setup_section():
    """Workspace setup section"""
    return [
        "",
        "RUN mkdir -p /workspace/logs && \\",
        "    chmod 777 /workspace"
    ]


def generate_copy_files_section(plugin_path):
    """Snakefile and scripts copy section"""
    lines = [
        "",
        "COPY Snakefile /workspace/Snakefile",
        "",
        "COPY scripts/ /scripts/",
        "",
        "WORKDIR /workspace"
    ]
    
    return lines


def generate_entrypoint_section(has_r=False):
    """Entrypoint script generation section - enhanced version"""
    lines = [
        "",
        "RUN echo '#!/bin/bash' > /entrypoint.sh && \\",
        "    echo 'export MAMBA_ROOT_PREFIX=/opt/micromamba' >> /entrypoint.sh && \\",
        "    echo 'export MAMBA_EXE=/usr/local/bin/micromamba' >> /entrypoint.sh && \\",
        "    echo 'export PATH=$MAMBA_ROOT_PREFIX/envs/plugin_env/bin:$PATH' >> /entrypoint.sh && \\"
    ]
    
    if has_r:
        lines.extend([
            "    echo 'export R_HOME=/opt/micromamba/envs/plugin_env/lib/R' >> /entrypoint.sh && \\",
            "    echo 'export R_LIBS_USER=/opt/micromamba/envs/plugin_env/lib/R/library' >> /entrypoint.sh && \\",
            "    echo 'export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE' >> /entrypoint.sh && \\"
        ])
    
    lines.extend([
        "    echo 'eval \"$($MAMBA_EXE shell activate -s bash -p $MAMBA_ROOT_PREFIX plugin_env)\" 2>/dev/null || true' >> /entrypoint.sh && \\",
        "    echo 'cd /workspace' >> /entrypoint.sh && \\",
        "    echo 'exec \"$@\"' >> /entrypoint.sh && \\",
        "    chmod +x /entrypoint.sh",
        "",
        "ENTRYPOINT [\"/entrypoint.sh\"]"
    ])
    
    return lines


def generate_healthcheck_section():
    """Health check section"""
    return [
        "",
        "HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \\",
        "    CMD test -f /opt/micromamba/envs/plugin_env/bin/python || exit 1"
    ]


def generate_cmd_section():
    """Default command section"""
    return [
        "",
        "CMD [\"/bin/bash\"]"
    ]


def generate_plugin_dockerfile(plugin_path: str, output_path: str, use_gpu: bool = True):
    """
    Analyzes plugin folder and generates multi-platform Dockerfile.
    
    Parameters:
        plugin_path (str): Plugin folder path
        output_path (str): Output Dockerfile path
        use_gpu (bool): GPU usage flag (default: True)
    """
    dependency_path = os.path.join(plugin_path, "dependency")
    
    # Check dependency files
    has_requirements = os.path.isfile(os.path.join(dependency_path, "requirements.txt"))
    has_environment = os.path.isfile(os.path.join(dependency_path, "environment.yml"))
    has_renv = os.path.isfile(os.path.join(dependency_path, "renv.lock"))
    
    # Check Python and R scripts
    has_python = False
    has_r = False
    scripts_path = os.path.join(plugin_path, "scripts")
    if os.path.exists(scripts_path):
        for file in os.listdir(scripts_path):
            if file.endswith(".py"):
                has_python = True
            if file.endswith(".R"):
                has_r = True
    
    # Extract versions
    python_version = "3.10"  # Default
    r_version = "4.4.2"  # Default
    
    if has_environment:
        python_version = extract_python_version_from_environment_yml(
            os.path.join(dependency_path, "environment.yml")
        )
    
    if has_renv:
        r_version = extract_r_version_from_renv_lock(
            os.path.join(dependency_path, "renv.lock")
        )
    
    # Generate Dockerfile
    dockerfile_lines = []
    
    # 1. Base image (multi-platform support)
    dockerfile_lines.extend(generate_base_image_section(use_gpu))
    
    # 2. Environment variables setup
    dockerfile_lines.extend(generate_env_setup_section())
    
    # 3. System packages installation
    dockerfile_lines.extend(generate_system_packages_section(use_gpu, has_r or has_renv))
    
    # 4. Dependency folder copy
    dockerfile_lines.extend(generate_dependency_copy_section())
    
    # 5. Micromamba installation (platform-specific)
    dockerfile_lines.extend(generate_micromamba_install_section(use_gpu))
    
    # 6. Environment variables configuration
    dockerfile_lines.extend(generate_env_variables_section(has_r or has_renv))
    
    # 7. Micromamba environment directory creation
    dockerfile_lines.extend(generate_micromamba_setup_section())
    
    # 8. Python environment creation
    dockerfile_lines.extend(generate_python_env_section(python_version))
    
    # 9. Snakemake and essential packages installation
    dockerfile_lines.extend(generate_snakemake_install_section())
    
    # 10. Python packages installation
    if has_requirements:
        dockerfile_lines.extend(generate_python_packages_section(has_requirements))
    
    # 11. R installation and packages (platform consideration)
    if has_r or has_renv:
        dockerfile_lines.extend(generate_r_installation_section(r_version, has_renv, use_gpu))
    
    # 12. Workspace setup
    dockerfile_lines.extend(generate_workspace_setup_section())
    
    # 13. Snakefile and scripts copy
    dockerfile_lines.extend(generate_copy_files_section(plugin_path))
    
    # 14. Entrypoint configuration
    dockerfile_lines.extend(generate_entrypoint_section(has_r or has_renv))
    
    # 15. Health check configuration
    dockerfile_lines.extend(generate_healthcheck_section())
    
    # 16. Default command configuration
    dockerfile_lines.extend(generate_cmd_section())
    
    # Save Dockerfile
    with open(output_path, "w") as f:
        f.write("\n".join(dockerfile_lines))
    
    # Output result logs
    print(f"[✓] Multi-platform Dockerfile generated for {plugin_path}")
    print(f"    - GPU: {use_gpu}")
    print(f"    - Platform support: {'linux/amd64' if use_gpu else 'linux/amd64, linux/arm64, linux/arm/v7'}")
    print(f"    - Python version: {python_version}")
    if has_r or has_renv:
        print(f"    - R version: {r_version}")
    print(f"    - has_requirements: {has_requirements}")
    print(f"    - has_renv: {has_renv}")
