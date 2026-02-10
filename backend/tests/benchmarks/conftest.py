"""
Benchmarks-specific conftest.

This file exists so that pytest can collect tests from this sub-package
without requiring the root conftest.py fixtures (which need httpx/TestClient).
No DB or HTTP fixtures are needed for function-level benchmarks.
"""
