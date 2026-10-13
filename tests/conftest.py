from pathlib import Path
import pytest
from cgcomp.fragments import load_fragments

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def fragments():
    return load_fragments(ROOT / 'single_nucleotide_fragments')
