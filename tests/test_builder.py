import anndata as ad
import numpy as np
import pytest
from scipy import sparse

from cellstate_c2s import DeviationSentenceBuilder

# Gene means over the 3 cells: A=2, B=2, C=1, D=2/3
X = np.array([[4, 0, 2, 0], [0, 3, 0, 1], [2, 3, 1, 1]], dtype=float)


def make_adata(dense=False):
    a = ad.AnnData(X.copy() if dense else sparse.csr_matrix(X))
    a.var_names = list("ABCD")
    a.obs_names = ["c0", "c1", "c2"]
    return a


@pytest.fixture(params=[False, True], ids=["sparse", "dense"])
def adata(request):
    return make_adata(dense=request.param)


def test_reference_mean(adata):
    b = DeviationSentenceBuilder().fit(adata)
    np.testing.assert_allclose(b.reference_mean_, [2, 2, 1, 2 / 3])
    assert b.n_reference_cells_ == 3


def test_masked(adata):
    names, s = DeviationSentenceBuilder().fit_transform(adata)
    assert names == ["c0", "c1", "c2"]
    assert s == ["A C", "B D", "B D A C"]  # last cell: A/C tie, stable order


def test_unmasked(adata):
    _, s = DeviationSentenceBuilder(mask_unexpressed=False).fit_transform(adata)
    assert s[0] == "A C D B"


def test_top_k(adata):
    _, s = DeviationSentenceBuilder(top_k=2, mask_unexpressed=False).fit_transform(adata)
    assert s[0] == "A C"


def test_chunking_matches(adata):
    _, full = DeviationSentenceBuilder().fit_transform(adata)
    _, chunked = DeviationSentenceBuilder(chunk_size=1).fit_transform(adata)
    assert full == chunked


def test_fit_mask(adata):
    # Reference from c0 only: A=4, B=0, C=2, D=0
    b = DeviationSentenceBuilder(mask_unexpressed=False).fit(adata, mask=[True, False, False])
    np.testing.assert_allclose(b.reference_mean_, [4, 0, 2, 0])
    _, s = b.transform(adata)
    assert s[1] == "B D C A"  # B=3, D=1, C=-2, A=-4


def test_fit_mask_errors(adata):
    with pytest.raises(ValueError):
        DeviationSentenceBuilder().fit(adata, mask=[False, False, False])
    with pytest.raises(ValueError):
        DeviationSentenceBuilder().fit(adata, mask=[True, False])


def test_transform_on_other_dataset_uses_shared_genes(adata):
    b = DeviationSentenceBuilder(mask_unexpressed=False).fit(adata)
    other = ad.AnnData(np.array([[5.0, 0.0, 9.0]]))
    other.var_names = ["C", "A", "Z"]  # different order, B and D missing, Z extra
    other.obs_names = ["x0"]
    with pytest.warns(UserWarning, match="2 of 4 reference genes"):
        names, s = b.transform(other)
    assert names == ["x0"]
    assert s == ["C A"]  # C: 5-1=4, A: 0-2=-2


def test_no_shared_genes(adata):
    b = DeviationSentenceBuilder().fit(adata)
    other = ad.AnnData(np.ones((1, 1)))
    other.var_names = ["Z"]
    with pytest.raises(ValueError, match="No genes shared"):
        b.transform(other)


def test_duplicate_var_names_rejected():
    a = ad.AnnData(np.ones((2, 2)))
    a.var_names = ["A", "A"]
    with pytest.raises(ValueError, match="unique"):
        DeviationSentenceBuilder().fit(a)


def test_transform_before_fit(adata):
    with pytest.raises(RuntimeError):
        DeviationSentenceBuilder().transform(adata)


def test_save_load_roundtrip(adata, tmp_path):
    b = DeviationSentenceBuilder(top_k=3, mask_unexpressed=False, chunk_size=7).fit(adata)
    b.save(tmp_path / "ref.npz")
    b2 = DeviationSentenceBuilder.load(tmp_path / "ref.npz")
    assert (b2.top_k, b2.mask_unexpressed, b2.chunk_size) == (3, False, 7)
    assert b.transform(adata) == b2.transform(adata)
