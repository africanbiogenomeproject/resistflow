"""
Tests for local data loading (data/local.py).
"""

import pytest
from resistflow.data.local import load_local_data


def _write_csv(tmp_path, content, name="data.csv"):
    path = tmp_path / name
    path.write_text(content)
    return path


def test_load_csv(tmp_path):
    path = _write_csv(tmp_path, "year,frequency,count,nobs\n2004,0.077,2,26\n2012,1.0,198,198\n")
    data = load_local_data(path)
    assert data["p0"] == pytest.approx(0.077)
    assert data["observed_years"] == [8.0]
    assert data["observed_frequencies"] == [1.0]
    assert data["n_chromosomes"] == 198
    assert data["p0_count"] == 2
    assert data["p0_nobs"] == 26
    assert data["baseline_year"] == 2004


def test_years_relative_to_baseline(tmp_path):
    path = _write_csv(tmp_path, "year,frequency\n2005,0.20\n2009,0.52\n2013,0.93\n")
    data = load_local_data(path)
    assert data["observed_years"] == [4.0, 8.0]


def test_unsorted_input(tmp_path):
    """Rows out of order are sorted; earliest year becomes the baseline."""
    path = _write_csv(tmp_path, "year,frequency\n2013,0.93\n2005,0.20\n2009,0.52\n")
    data = load_local_data(path)
    assert data["baseline_year"] == 2005
    assert data["p0"] == pytest.approx(0.20)
    assert data["observed_years"] == [4.0, 8.0]


def test_column_aliases_and_case(tmp_path):
    """Alternative column names and capitalisation are accepted."""
    path = _write_csv(tmp_path, "Sampling_Year,AF\n2004,0.077\n2012,1.0\n")
    data = load_local_data(path)
    assert data["p0"] == pytest.approx(0.077)
    assert data["observed_years"] == [8.0]


def test_optional_columns_absent(tmp_path):
    path = _write_csv(tmp_path, "year,frequency\n2004,0.077\n2012,1.0\n")
    data = load_local_data(path)
    assert data["n_chromosomes"] is None
    assert data["p0_count"] is None
    assert data["p0_nobs"] is None


def test_missing_required_column_raises(tmp_path):
    path = _write_csv(tmp_path, "year,value\n2004,0.077\n2012,1.0\n")
    with pytest.raises(ValueError) as excinfo:
        load_local_data(path)
    assert "frequency" in str(excinfo.value)


def test_single_observation_raises(tmp_path):
    path = _write_csv(tmp_path, "year,frequency\n2004,0.077\n")
    with pytest.raises(ValueError):
        load_local_data(path)


def test_unsupported_file_type_raises(tmp_path):
    path = _write_csv(tmp_path, "year,frequency\n2004,0.077\n", name="data.json")
    with pytest.raises(ValueError):
        load_local_data(path)


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_local_data(tmp_path / "does_not_exist.csv")