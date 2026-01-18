
import subprocess
import pytest

def test_promoter_retrieve_help():
    result = subprocess.run(["promoter-retrieve", "--help"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "Usage" in result.stdout

def test_get_genome_statistic_help():
    result = subprocess.run(["get-genome-statistic", "--help"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "Usage" in result.stdout

def test_proteins_retrieve_help():
    result = subprocess.run(["proteins-retrieve", "--help"], capture_output=True, text=True)
    assert result.returncode == 0
    assert "Usage" in result.stdout
