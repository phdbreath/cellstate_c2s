"""
Glue for Cell2Sentence arrow datasets.

Works on Hugging Face datasets.Dataset, such as the one returned by
cell2sentence.CSData.adata_to_arrow. This module does not import
cell2sentence itself.
"""


def replace_sentences(arrow_ds, cell_names, sentences,
                      field="cell_sentence", name_field="cell_name"):
    """Swap the sentence column of an arrow dataset, matching rows by cell name.

    Parameters
    ----------
    arrow_ds : datasets.Dataset
        Dataset with (field) and (name_field) columns.
    cell_names, sentences : sequence of str
        Output of (DeviationSentenceBuilder.transform).
    field : str
        Column holding the cell sentence.
    name_field : str
        Column holding the cell name, matching (adata.obs_names).

    Returns
    -------
    datasets.Dataset
        A new dataset; the input is not modified.
    """
    if len(cell_names) != len(sentences):
        raise ValueError("cell_names and sentences must have the same length.")
    lookup = dict(zip(cell_names, sentences))
    if len(lookup) != len(cell_names):
        raise ValueError("cell_names contains duplicates.")
    for col in (field, name_field):
        if col not in arrow_ds.column_names:
            raise KeyError(f"Column {col!r} not in dataset; columns are {arrow_ds.column_names}.")
    missing = [n for n in arrow_ds[name_field] if n not in lookup]
    if missing:
        raise KeyError(
            f"{len(missing)} cells in the arrow dataset have no sentence, e.g. {missing[:3]}"
        )
    return arrow_ds.map(lambda ex: {field: lookup[ex[name_field]]})
