"""Shared pytest fixtures: locate repo data directories and load fixture files."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = REPO_ROOT / "examples"
VALIDATION_DIR = REPO_ROOT / "validation"


@pytest.fixture(scope="session")
def examples_dir() -> Path:
    return EXAMPLES_DIR


@pytest.fixture(scope="session")
def validation_dir() -> Path:
    return VALIDATION_DIR


def load_example(name: str) -> str:
    return (EXAMPLES_DIR / name).read_text(encoding="utf-8")


def load_validation(name: str) -> dict:
    return json.loads((VALIDATION_DIR / name).read_text(encoding="utf-8"))
