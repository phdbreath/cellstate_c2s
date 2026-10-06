import pytest

datasets = pytest.importorskip("datasets")

from cellstate_c2s import replace_sentences  # noqa: E402


@pytest.fixture
def arrow_ds():
    return datasets.Dataset.from_dict({
        "cell_name": ["c1", "c0"],
        "cell_sentence": ["old1", "old0"],
        "annotation_cell_states": ["AT2_2", "AT2_1"],
    })


def test_replace_matches_by_name(arrow_ds):
    out = replace_sentences(arrow_ds, ["c0", "c1"], ["A C", "B D"])
    assert out["cell_sentence"] == ["B D", "A C"]
    assert out["annotation_cell_states"] == ["AT2_2", "AT2_1"]
    assert arrow_ds["cell_sentence"] == ["old1", "old0"]  # input untouched


def test_missing_cell_raises(arrow_ds):
    with pytest.raises(KeyError, match="no sentence"):
        replace_sentences(arrow_ds, ["c0"], ["A C"])


def test_duplicate_names_raise(arrow_ds):
    with pytest.raises(ValueError, match="duplicates"):
        replace_sentences(arrow_ds, ["c0", "c0", "c1"], ["a", "b", "c"])


def test_missing_column_raises(arrow_ds):
    with pytest.raises(KeyError, match="not in dataset"):
        replace_sentences(arrow_ds, ["c0", "c1"], ["a", "b"], name_field="barcode")
