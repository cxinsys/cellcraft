"""
Plugin Docker build helpers.

Extracted from ``plugin/router.py`` (PR-5). This module concentrates the
Docker/Dockerfile-facing logic that the build endpoints used to run inline:

- ``prepare_build_context()``: ensure the scripts folder exists (with a dummy
  ``.gitkeep`` when empty) and generate the plugin's Dockerfile.
- ``dispatch_build_task()``: kick off the asynchronous Celery build task.

The heavy lifting still lives in ``plugin/utils.py`` (``generate_plugin_dockerfile``)
and ``app.worker.tasks`` (``build_plugin_task``); ``utils.py`` internals are not
touched here (that split is PR-9). This module only owns the orchestration that
previously sat in the router.
"""
import os

from app.plugin import utils as plugin_utils
from app.worker.tasks import build_plugin_task


def prepare_build_context(*, plugin_folder: str, script_folder: str, use_gpu: bool) -> str:
    """Ensure the scripts folder is build-ready and generate the Dockerfile.

    Mirrors the inline logic from ``build_plugin_docker``:
    - create the scripts folder if missing
    - drop a ``.gitkeep`` dummy file when the scripts folder is empty (Docker
      build needs a non-empty directory)
    - generate the Dockerfile at ``<plugin_folder>/Dockerfile``

    Returns the generated Dockerfile path.
    """
    # scripts 폴더 존재 확인 및 더미 파일 생성 (Docker 빌드를 위해 필요)
    print(f"Checking scripts folder before Docker build: {script_folder}")
    print(f"Scripts folder exists: {os.path.exists(script_folder)}")

    if not os.path.exists(script_folder):
        os.makedirs(script_folder)
        print(f"Created empty scripts folder at {script_folder}")

    # scripts 폴더가 비어있다면 더미 파일 생성
    scripts_contents = os.listdir(script_folder) if os.path.exists(script_folder) else []
    print(f"Scripts folder contents: {scripts_contents}")

    if not scripts_contents:
        dummy_file_path = os.path.join(script_folder, ".gitkeep")
        with open(dummy_file_path, 'w') as f:
            f.write("# This file ensures the scripts directory is not empty\n")
        print(f"Created dummy file at {dummy_file_path}")
        print(f"Scripts folder contents after dummy file: {os.listdir(script_folder)}")

    # Dockerfile 생성
    dockerfile_path = os.path.join(plugin_folder, "Dockerfile")
    plugin_utils.generate_plugin_dockerfile(plugin_folder, dockerfile_path, use_gpu=use_gpu)
    print(f"Generated Dockerfile at: {dockerfile_path} (GPU: {use_gpu})")

    return dockerfile_path


def dispatch_build_task(*, plugin_name: str, user_id: int):
    """Dispatch the asynchronous plugin Docker build Celery task.

    Wire format (kwargs keys) is preserved exactly as the router used it.
    Returns the Celery ``AsyncResult``.
    """
    return build_plugin_task.apply_async(
        args=[],
        kwargs={
            'plugin_name': plugin_name,
            'user_id': user_id,
            'workflow_id': None,
            'algorithm_id': None,
            'task_type': "plugin_build"
        },
        ignore_result=False
    )
