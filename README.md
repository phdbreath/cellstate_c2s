# cellstate_c2s

Custom cell sentences for [Cell2Sentence](https://github.com/vandijklab/cell2sentence) fine-tuning and prediction.

Cell2Sentence normally writes each cell as its genes ordered by expression. This package replaces that with genes ordered by **deviation from a reference mean** (for example, the mean AT2 profile), so the sentence describes how a cell differs from the reference state rather than which genes are simply most abundant. You build the arrow dataset with Cell2Sentence as usual, then swap in the new sentences before creating `CSData`.

## Install

```bash
pip install cellstate_c2s            # core: numpy, scipy, pandas, anndata
pip install "cellstate_c2s[c2s]"     # plus cell2sentence and datasets
```

From a local clone:

```bash
pip install -e ".[dev]"
pytest
```

## Usage

```python
import cell2sentence as cs
from cellstate_c2s import DeviationSentenceBuilder, replace_sentences

# 1. Build the arrow dataset with Cell2Sentence as usual
arrow_ds, vocabulary = cs.CSData.adata_to_arrow(
    adata=adata_train, random_state=SEED, sentence_delimiter=" ",
    label_col_names=["annotation_cell_states", "organism", "tissue"],
)

# 2. Fit the reference and build deviation sentences
builder = DeviationSentenceBuilder(top_k=200)
builder.fit(adata_train, mask=adata_train.obs["cell_type"] == "AT2")  # optional mask
names, sentences = builder.transform(adata_train)
builder.save("reference.npz")

# 3. Swap the sentences in, matched by cell name, then continue with Cell2Sentence
arrow_ds = replace_sentences(arrow_ds, names, sentences)
csdata = cs.CSData.csdata_from_arrow(
    arrow_dataset=arrow_ds, vocabulary=vocabulary,
    save_dir=save_dir, save_name=save_name, dataset_backend="arrow",
)
```

**Held-out test data (same platform):** reuse the reference fitted on training data. Do not fit on the test set.

```python
builder = DeviationSentenceBuilder.load("reference.npz")
names, sentences = builder.transform(adata_test)
```

**Another platform (e.g. Xenium):** fit a separate reference on that dataset, then do the Cell2Sentence vocabulary swap yourself.

```python
names_x, sents_x = DeviationSentenceBuilder(top_k=200).fit_transform(adata_xenium)
arrow_ds_xenium = replace_sentences(arrow_ds_xenium, names_x, sents_x)
```

## How sentences are built

For each cell, `deviation = expression - reference_mean`. Genes are sorted by deviation, highest first, and the first `top_k` are kept. Ties keep the order of `var_names` (stable sort), so the same input always gives the same sentence.

| Option | Default | Effect |
|---|---|---|
| `top_k` | 200 | Maximum genes per sentence. Capped at the number of shared genes. |
| `mask_unexpressed` | `True` | Drops genes whose value in that cell is exactly 0, so sentences can be shorter than `top_k`. Set `False` for a pure ranking of all genes. |
| `chunk_size` | 2048 | Cells densified at a time. Lower it if memory is tight. |
| `fit(..., mask=)` | all cells | Boolean array selecting the cells that define the reference. |

`transform` uses only the genes shared by the reference and the AnnData passed in, warns about missing reference genes, and ignores extra genes.

## Before you use this

The code cannot check most of these, so they are on you. The first four are the most likely to cause silent errors.

1. **`adata.X` is normalized, the same way in every dataset.** Deviation is in expression units, so raw counts in one dataset and log-normalized values in another give meaningless ranks. The package never reads `layers` or `.raw`.
2. **Zero means "not expressed."** With `mask_unexpressed=True`, any gene with a value of exactly 0 in a cell is dropped.
3. **The cells passed to `fit` define the reference.** Pass only the cells you want (e.g. AT2), or use `mask=`.
4. **Fit on training data only for evaluation.** `fit_transform` on a test set computes the reference from test cells.
5. **Gene names are unique and use one naming system** (symbols or Ensembl IDs, not a mix). Duplicates raise an error; mismatched naming just shrinks the shared gene set.
6. **Cell names (`obs_names`) are unique and match the arrow dataset's `cell_name` column.** Field names are parameters of `replace_sentences` if your Cell2Sentence version differs.
7. **The datasets are already subset to the genes you want.** The package does not select genes.

## What is checked

| Checked in code | Not checked |
|---|---|
| `fit` called before `transform` | Normalization |
| Unique gene names | Gene naming system |
| Shared genes exist; warns on missing ones | Whether the reference cells are the right ones |
| Every arrow row has a sentence; no duplicate cell names | Train/test separation |
| Required arrow columns exist | Meaning of zeros |

## License

MIT
