"""Deviation-from-reference cell sentences."""

import warnings

import numpy as np
import pandas as pd
from scipy import sparse


class DeviationSentenceBuilder:
    """
    Build cell sentences by ranking genes on signed deviation from a reference mean.

    For each cell, (deviation = expression - reference_mean). Genes are sorted
    by deviation (highest first, stable on ties) and the first (top_k) are
    kept. With (mask_unexpressed=True), genes whose value in that cell is
    exactly 0 are dropped, so sentences can be shorter than (top_k).

    Parameters
    ----------
    top_k : int
        Maximum number of genes per sentence.
    mask_unexpressed : bool
        Drop genes with zero expression in the cell.
    chunk_size : int
        Number of cells densified at a time in ``transform``.

    Notes
    -----
    adata.X must hold normalized expression, on the same scale for every
    dataset you convert. The reference mean is taken over every cell passed to
    fit (or the subset selected by mask).
    """

    def __init__(self, top_k=200, mask_unexpressed=True, chunk_size=2048):
        if top_k < 1:
            raise ValueError("top_k must be at least 1.")
        if chunk_size < 1:
            raise ValueError("chunk_size must be at least 1.")
        self.top_k = top_k
        self.mask_unexpressed = mask_unexpressed
        self.chunk_size = chunk_size
        self.reference_mean_ = None
        self.gene_names_ = None
        self.n_reference_cells_ = None

    def fit(self, adata, mask=None):
        """
        Compute the reference mean.

        Parameters
        ----------
        adata : anndata.AnnData
            Cells x genes, normalized.
        mask : array-like of bool, optional
            Selects the cells that define the reference (for example AT2 cells
            only). Defaults to all cells.
        """
        _check_unique_var_names(adata)
        X = adata.X
        if mask is not None:
            mask = np.asarray(mask, dtype=bool)
            if mask.shape != (adata.n_obs,):
                raise ValueError(f"mask must have shape ({adata.n_obs},), got {mask.shape}.")
            if not mask.any():
                raise ValueError("mask selects no cells.")
            X = X[mask]
        if X.shape[0] == 0:
            raise ValueError("Cannot fit a reference on zero cells.")
        self.reference_mean_ = np.asarray(X.mean(axis=0)).ravel().astype(np.float64)
        self.gene_names_ = np.asarray(adata.var_names).astype(str)
        self.n_reference_cells_ = int(X.shape[0])
        return self

    def transform(self, adata):
        """
        Return (cell_names, sentences) for every cell in adata.

        Only genes present in both the reference and adata are used. A
        warning lists how many reference genes are missing.
        """
        self._check_fitted()
        _check_unique_var_names(adata)
        idx = pd.Index(adata.var_names).get_indexer(self.gene_names_)
        keep = idx >= 0
        if not keep.any():
            raise ValueError("No genes shared between the reference and this AnnData.")
        if not keep.all():
            warnings.warn(
                f"{int((~keep).sum())} of {len(keep)} reference genes are missing from "
                "this AnnData and will be ignored.",
                stacklevel=2,
            )
        col_idx = idx[keep]
        genes = self.gene_names_[keep]
        ref = self.reference_mean_[keep]
        k = min(self.top_k, len(genes))

        X = adata.X
        sentences = []
        for start in range(0, X.shape[0], self.chunk_size):
            chunk = X[start:start + self.chunk_size][:, col_idx]
            chunk = chunk.toarray() if sparse.issparse(chunk) else np.asarray(chunk)
            chunk = chunk.astype(np.float64, copy=False)
            dev = chunk - ref
            if self.mask_unexpressed:
                dev = np.where(chunk > 0, dev, -np.inf)
            order = np.argsort(-dev, axis=1, kind="stable")[:, :k]
            for row in range(order.shape[0]):
                o = order[row]
                if self.mask_unexpressed:
                    o = o[np.isfinite(dev[row, o])]
                sentences.append(" ".join(genes[o]))
        return [str(n) for n in adata.obs_names], sentences

    def fit_transform(self, adata, mask=None):
        """Fit the reference on ``adata`` and transform the same cells."""
        return self.fit(adata, mask=mask).transform(adata)

    def save(self, path):
        """Save the reference and settings to a ``.npz`` file."""
        self._check_fitted()
        np.savez(
            path,
            reference_mean=self.reference_mean_,
            gene_names=self.gene_names_,
            n_reference_cells=self.n_reference_cells_,
            top_k=self.top_k,
            mask_unexpressed=self.mask_unexpressed,
            chunk_size=self.chunk_size,
        )

    @classmethod
    def load(cls, path):
        """Load a builder saved with :meth:save."""
        with np.load(path, allow_pickle=False) as d:
            obj = cls(
                top_k=int(d["top_k"]),
                mask_unexpressed=bool(d["mask_unexpressed"]),
                chunk_size=int(d["chunk_size"]),
            )
            obj.reference_mean_ = d["reference_mean"]
            obj.gene_names_ = d["gene_names"].astype(str)
            obj.n_reference_cells_ = int(d["n_reference_cells"])
        return obj

    def _check_fitted(self):
        if self.reference_mean_ is None:
            raise RuntimeError("Call fit() first.")


def _check_unique_var_names(adata):
    if not adata.var_names.is_unique:
        raise ValueError("var_names must be unique; call adata.var_names_make_unique().")
