import numpy as np
import pytest
from cgcomp.parsers.oxrna import parse, read_topology
from conftest import ROOT

DAT = ROOT / 'CG_models/exp_oxRNA.dat'
TOP = ROOT / 'examples/auxiliary/exp_oxRNA.top'


def test_actual_oxrna_native_configuration_and_topology():
    rows = np.loadtxt(DAT, skiprows=3)
    parsed = parse(DAT, TOP)
    assert len(parsed) == 102
    # Compare sequence with the authentic companion topology, independently.
    assert ''.join(r.base for r in parsed) == TOP.read_text().splitlines()[1].split()[0]
    for r, row in zip(parsed, rows):
        assert np.array_equal(r.sites['COM'], row[:3])
        assert np.array_equal(r.orientation[0], row[3:6])
        assert np.array_equal(r.orientation[2], row[6:9])
        assert np.allclose(r.orientation @ r.orientation.T, np.eye(3), atol=1e-12)
        assert np.linalg.det(r.orientation) == pytest.approx(1)
    assert parsed[0].neighbors == (1, -1)
    assert parsed[-1].neighbors == (-1, 100)


def test_legacy_topology_equivalent_connectivity(tmp_path):
    records = read_topology(TOP)
    old = tmp_path / 'legacy.top'
    old.write_text('102 1\n' + ''.join(f'{s} {b} {n3} {n5}\n' for s, b, n3, n5 in records))
    assert read_topology(old) == records
    assert len(parse(DAT, old)) == 102


def test_oxrna_missing_topology():
    with pytest.raises(ValueError, match='matching native topology'):
        parse(DAT, None)


def test_invalid_oxrna_orientation(tmp_path):
    lines = DAT.read_text().splitlines()
    fields = lines[3].split()
    fields[3:6] = ['0', '0', '0']
    lines[3] = ' '.join(fields)
    invalid = tmp_path / 'invalid.dat'
    invalid.write_text('\n'.join(lines))
    with pytest.raises(ValueError, match='orthonormal'):
        parse(invalid, TOP)
