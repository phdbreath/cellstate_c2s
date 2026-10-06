"""End-to-end example: deviation sentences for Cell2Sentence fine-tuning and prediction.

Edit the paths and column names, then run cell by cell in a notebook or as a script.
Requires: pip install "cell-state-c2s[c2s]" scanpy

For fine-tuning and normalization related information, please check the cell2sentence documentation: https://vandijklab-cell2sentence.readthedocs.io/en/latest/index.html
"""

import cell2sentence as cs
import scanpy as sc

from cellstate_c2s import DeviationSentenceBuilder, replace_sentences

SEED = 1234
LABEL_COLS = ["annotation_cell_states", "organism", "tissue"]

# Fine-tuning dataset 
adata_train = sc.read_h5ad("fine_tune.h5ad")  # normalized, shared genes only
adata_train.obs["organism"] = "human"
adata_train.obs["tissue"] = "lung"
for col in LABEL_COLS:
    adata_train.obs[col] = adata_train.obs[col].astype(str)

arrow_train, vocabulary = cs.CSData.adata_to_arrow(
    adata=adata_train, random_state=SEED, sentence_delimiter=" ",
    label_col_names=LABEL_COLS,
)

builder = DeviationSentenceBuilder(top_k=200).fit(adata_train)
builder.save("reference_fine_tune.npz")
names, sentences = builder.transform(adata_train)
arrow_train = replace_sentences(arrow_train, names, sentences)
print(arrow_train[0])

csdata_train = cs.CSData.csdata_from_arrow(
    arrow_dataset=arrow_train, vocabulary=vocabulary,
    save_dir="c2s_data", save_name="fine_tune_train", dataset_backend="arrow",
)

# Fine-tune a model on the training data (check cell2sentence docs to understand how)
# Prediction on the target dataset 
adata_target = sc.read_h5ad("target.h5ad")  # same normalization as training
adata_target.obs["organism"] = "human"
adata_target.obs["tissue"] = "lung"
for col in LABEL_COLS:
    adata_target.obs[col] = adata_target.obs[col].astype(str)

arrow_target, vocabulary_target = cs.CSData.adata_to_arrow(
    adata=adata_target, random_state=SEED, sentence_delimiter=" ",
    label_col_names=LABEL_COLS,
)

# target gets its own reference; then swap the vocabulary as your C2S workflow requires.
names_x, sents_x = DeviationSentenceBuilder(top_k=200).fit_transform(adata_target)
arrow_target = replace_sentences(arrow_target, names_x, sents_x)

csdata_target = cs.CSData.csdata_from_arrow(
    arrow_dataset=arrow_target, vocabulary=vocabulary_target,
    save_dir="c2s_data", save_name="target", dataset_backend="arrow",
)

# Predict the target dataset using the fine-tuned model
