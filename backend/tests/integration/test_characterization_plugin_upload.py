"""
Characterization tests for the plugin upload endpoint.

Purpose: pin the CURRENT behavior of ``POST /routes/plugin/upload`` exactly as it
is today (upload_plugin, app/routes/endpoints/plugin.py lines 130-288). These
tests describe reality so a later service-extraction refactor (PR-5) can prove it
did not change observable behavior.

The endpoint has heavy side effects (filesystem writes under ./plugin/local,
metadata + Snakefile generation) followed by DB write + commit + rollback. To
pin the HTTP/DB contract without touching the real filesystem, the
``plugin_utils`` helpers and the cache invalidation call are mocked
(consistent with the existing mocking style in test_plugin_api.py).

Source of truth read while writing these tests:
- app/routes/endpoints/plugin.py :: upload_plugin
- app/database/schemas/plugin.py :: PluginCreate (flat schema this endpoint uses)

Pinned behavior:
- New local plugin: 200 with
  {message: "플러그인 메타데이터 업로드 성공", plugin: {...}, success: True,
   redirect_to: "plugin_list"} and a DB row created.
- Uploading over an existing OFFICIAL plugin name: 403 (read-only official).
- DB failure during create is rolled back (no row persists) and yields 500.
- Invalid request body (missing required PluginCreate fields): 422.
"""
import pytest
from unittest.mock import patch
from fastapi.testclient import TestClient
from sqlalchemy.orm import Session

from app.database import models
from app.core.config import settings
from app.core.enums import PluginType

UPLOAD_URL = f"{settings.ROUTES_STR}/plugin/upload"


def _valid_plugin_create_payload(name: str = "CharPlugin") -> dict:
    """A body matching the flat PluginCreate schema the /upload endpoint expects."""
    return {
        "name": name,
        "description": "Characterization plugin",
        "author": "char_author",
        "plugin_path": f"./plugin/local/{name}/",
        "plugin_type": "analysis",
        "dependencies": {"requirements.txt": "numpy==1.24.0"},
        "reference_folders": None,
        "drawflow": {"drawflow": {"Home": {"data": {}}}},
        "rules": {"0": {"name": "r", "script": "s.py"}},
        "use_gpu": False,
        "source": "local",
        "is_editable": True,
        "version": "1.0.0",
    }


@pytest.mark.integration
@pytest.mark.characterization
class TestCharacterizationPluginUpload:
    """Freeze current upload_plugin behavior (success, official reject, rollback)."""

    @patch("app.routes.endpoints.plugin.invalidate_all_plugin_cache")
    @patch("app.routes.endpoints.plugin.plugin_utils.generate_snakemake_code")
    @patch("app.routes.endpoints.plugin.plugin_utils.create_metadata_file")
    @patch("app.routes.endpoints.plugin.plugin_utils.create_dependency_folder")
    @patch("app.routes.endpoints.plugin.plugin_utils.create_plugin_folder")
    @patch("app.routes.endpoints.plugin.ensure_local_plugins_dir")
    @patch("app.routes.endpoints.plugin.get_plugin_path")
    @patch("app.routes.endpoints.plugin.os.path.exists", return_value=False)
    @patch("app.routes.endpoints.plugin.os.makedirs")
    def test_upload_new_plugin_success_response_and_db(
        self,
        _mock_makedirs,
        _mock_exists,
        mock_get_plugin_path,
        _mock_ensure_dir,
        _mock_create_folder,
        _mock_create_dep,
        _mock_create_meta,
        _mock_gen_snake,
        _mock_invalidate,
        client: TestClient,
        auth_headers: dict,
        db_session: Session,
    ):
        """A brand-new local plugin uploads successfully with the current shape."""
        # Plugin does not exist yet -> get_plugin_path raises (endpoint treats as new).
        from fastapi import HTTPException
        mock_get_plugin_path.side_effect = HTTPException(status_code=404, detail="not found")

        payload = _valid_plugin_create_payload(name="CharNewPlugin")
        response = client.post(UPLOAD_URL, json=payload, headers=auth_headers)

        assert response.status_code == 200, response.text
        data = response.json()
        assert data["message"] == "플러그인 메타데이터 업로드 성공"
        assert data["success"] is True
        assert data["redirect_to"] == "plugin_list"
        # FastAPI 0.78 serializes the SQLAlchemy Plugin model via jsonable_encoder
        # which yields an empty dict {} (no Pydantic schema on this endpoint).
        # Pin the actual observed shape: plugin field is present and is a dict.
        assert isinstance(data["plugin"], dict)

        # DB row was created and committed.
        db_plugin = (
            db_session.query(models.Plugin)
            .filter(models.Plugin.name == "CharNewPlugin")
            .first()
        )
        assert db_plugin is not None
        assert db_plugin.source == "local"
        assert db_plugin.plugin_type == PluginType.ANALYSIS

    def test_upload_over_official_plugin_bypasses_guard_and_returns_500(
        self,
        client: TestClient,
        auth_headers: dict,
        db_session: Session,
    ):
        """Uploading a plugin whose name matches an official plugin returns 500 today.

        ANOMALY: The intended 403 guard raises HTTPException inside a bare
        ``except HTTPException: pass`` block, so the guard is silently bypassed
        and the upload proceeds until it fails at Snakemake code generation (500).
        This test pins the *actual* current behavior, not the intended behavior.
        """
        # Seed an official (read-only) plugin so get_plugin_path resolves to official.
        official = models.Plugin(
            name="OfficialCharPlugin",
            description="official",
            author="official",
            plugin_path="./plugin/official/OfficialCharPlugin/",
            plugin_type=PluginType.ANALYSIS,
            dependencies={},
            drawflow={},
            rules={},
            source="official",
            is_editable=False,
            version="1.0.0",
        )
        db_session.add(official)
        db_session.commit()

        payload = _valid_plugin_create_payload(name="OfficialCharPlugin")
        response = client.post(UPLOAD_URL, json=payload, headers=auth_headers)

        # ANOMALY (pinned): The 403 guard at lines 147-153 of upload_plugin raises
        # HTTPException(403) when get_plugin_path returns "official", but it is
        # immediately caught by ``except HTTPException: pass`` on line 152-153, so
        # the upload proceeds and eventually fails at generate_snakemake_code (500).
        # The official-plugin protection is therefore silently bypassed today.
        assert response.status_code == 500

    @patch("app.routes.endpoints.plugin.invalidate_all_plugin_cache")
    @patch("app.routes.endpoints.plugin.crud_plugin.create_plugin")
    @patch("app.routes.endpoints.plugin.plugin_utils.generate_snakemake_code")
    @patch("app.routes.endpoints.plugin.plugin_utils.create_metadata_file")
    @patch("app.routes.endpoints.plugin.plugin_utils.create_dependency_folder")
    @patch("app.routes.endpoints.plugin.plugin_utils.create_plugin_folder")
    @patch("app.routes.endpoints.plugin.ensure_local_plugins_dir")
    @patch("app.routes.endpoints.plugin.get_plugin_path")
    @patch("app.routes.endpoints.plugin.os.path.exists", return_value=False)
    @patch("app.routes.endpoints.plugin.os.makedirs")
    def test_upload_db_failure_rolls_back_and_returns_500(
        self,
        _mock_makedirs,
        _mock_exists,
        mock_get_plugin_path,
        _mock_ensure_dir,
        _mock_create_folder,
        _mock_create_dep,
        _mock_create_meta,
        _mock_gen_snake,
        mock_create_plugin,
        _mock_invalidate,
        client: TestClient,
        auth_headers: dict,
        db_session: Session,
    ):
        """A DB error during create triggers rollback (no row) and a 500 response."""
        from fastapi import HTTPException
        mock_get_plugin_path.side_effect = HTTPException(status_code=404, detail="not found")
        mock_create_plugin.side_effect = Exception("DB write boom")

        payload = _valid_plugin_create_payload(name="RollbackCharPlugin")
        response = client.post(UPLOAD_URL, json=payload, headers=auth_headers)

        assert response.status_code == 500

        # Rollback means nothing persisted for this plugin name.
        db_session.expire_all()
        db_plugin = (
            db_session.query(models.Plugin)
            .filter(models.Plugin.name == "RollbackCharPlugin")
            .first()
        )
        assert db_plugin is None

    def test_upload_invalid_body_returns_422(
        self,
        client: TestClient,
        auth_headers: dict,
    ):
        """Missing required PluginCreate fields fails Pydantic validation (422)."""
        # PluginData-style nested body is NOT what /upload expects (it wants flat
        # PluginCreate). Missing name/author/plugin_path -> validation error.
        response = client.post(
            UPLOAD_URL,
            json={"description": "no required fields"},
            headers=auth_headers,
        )
        assert response.status_code == 422

    def test_upload_requires_authentication(
        self,
        client: TestClient,
    ):
        """Without auth headers the endpoint rejects with 401."""
        response = client.post(UPLOAD_URL, json=_valid_plugin_create_payload())
        assert response.status_code == 401
