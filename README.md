# cellstate_c2s

Cell2Sentence 'sentences' are redefined for fine tuning and prediction to find cell states accurately.

Credits/Disclaimer: The original package is [Cell2Sentence](https://github.com/vandijklab/cell2sentence) . Our approach is modification to highlight cell state differences.


Cell2Sentence creates a cell sentence for every cell. The order of the words is based on the expression level. It orders genes from highest expressing to the lowest expressing. While this method works well for cell types, it struggles for cell states. The reason is intuitive. Cell states share the same high expressing genes and trying to differentiate the cell states becomes a difficult task because the cell state defining features are not placed in the early positions of the sentence. We propose reframing the solution. For every cell type independently, we calculated the mean of every gene in that cell type. Then, we look at how the expression is deviated from the mean cell type signature. This provides a modified cell sentence where important markers of the cell states are being placed in the beginning of the sentence. This increases the accuracy of the predictions.


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
| `top_k` | 200 | Maximum genes per sentence. |
| `mask_unexpressed` | `True` | Drops genes whose value in that cell is exactly 0, so sentences can be shorter than `top_k`. Set `False` for a pure ranking of all genes. The original sentence drops zero expressing cells. In our experiments we tested both versions for the modified model. We also recommend dropping them for cell states. |
| `chunk_size` | 2048 | Cells densified at a time. (Most of the time it is not a bottleneck) |
| `fit(..., mask=)` | all cells | Boolean array selecting the cells that define the reference. |

`transform` uses only the genes shared by the reference and the AnnData passed in, warns about missing reference genes, and ignores extra genes.

## Before you use this

This package modifies cell sentences. You need to implement recommended steps on the cell2sentence documentation. After you create the sentences, you can use this package to modify the sentences. There are some assumptions that this package makes.

1. **`adata.X` is normalized, the same way in every dataset.** Cell2Sentence requires the normalized values, this assumes it is done and you made sure X holds the normalized values.
2. **Zero means "not expressed."** With `mask_unexpressed=True`, any gene with a value of exactly 0 in a cell is dropped [Recommended to drop].
3. **The cells passed to `fit` define the reference.** Pass only the cells you want (e.g. AT2), or use `mask=`.
4. **Fit on training data only for evaluation.** `fit_transform` on a test set computes the reference from test cells.
5. **Gene names are unique and use one naming system** Both the finetuning and prediction datasets have the same gene names and you are using the common genes on both.
6. **Cell names (`obs_names`) are unique and match the arrow dataset's `cell_name` column.** Field names are parameters of `replace_sentences` if your Cell2Sentence version differs.
7. **The datasets are already subset to the genes you want.** You must subset the cell type of interest before finetuning, and also for prediction.

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
