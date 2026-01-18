
import pytest
from genomemanagement.seq_manage import FastaManager

@pytest.fixture
def mock_fasta(tmp_path):
    d = tmp_path / "subdir"
    d.mkdir()
    p = d / "test.fa"
    p.write_text(">seq1\nATGCATGC\n>seq2\nGGCCGGCC\n")
    return str(p)

def test_fasta_manager_init(mock_fasta):
    fm = FastaManager(mock_fasta)
    assert "seq1" in fm.chromosome_seq
    assert "seq2" in fm.chromosome_seq
    assert fm.chromosome_seq["seq1"] == "ATGCATGC"

def test_get_gc_content(mock_fasta):
    fm = FastaManager(mock_fasta)
    gc = fm.get_gc_content("ATGC")
    assert gc == 50.0
    gc = fm.get_gc_content("GGCC")
    assert gc == 100.0

def test_complementary(mock_fasta):
    fm = FastaManager(mock_fasta)
    seq = "ATGC"
    comp = fm.complementary(seq)
    assert comp == "TACG"
