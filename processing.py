"""Gene expression preprocessing — mirrors WEKA Normalize + Ranker (F-score)."""

from __future__ import annotations

import io

import numpy as np
import pandas as pd
from sklearn.feature_selection import f_classif


def _gene_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in df.columns if col.startswith("Gene_")]


def _detect_label_column(df: pd.DataFrame) -> str:
    non_gene = [col for col in df.columns if not col.startswith("Gene_")]
    if not non_gene:
        raise ValueError(
            "No label column found. Expected a column like 'Diagnosis' alongside Gene_* columns."
        )
    if "Diagnosis" in non_gene:
        return "Diagnosis"
    return non_gene[-1]


def load_dataset(file_bytes: bytes, filename: str) -> pd.DataFrame:
    """Load uploaded CSV or ARFF file."""
    name = filename.lower()
    if name.endswith(".arff"):
        from scipy.io import arff

        data, _meta = arff.loadarff(io.BytesIO(file_bytes))
        df = pd.DataFrame(data)
        for col in df.select_dtypes(include=["object"]).columns:
            df[col] = df[col].astype(str).str.replace("b'", "", regex=False).str.replace("'", "", regex=False)
        return df

    return pd.read_csv(io.BytesIO(file_bytes))


def clean_dataset(df: pd.DataFrame, label_col: str | None = None) -> tuple[pd.DataFrame, dict]:
    """
    Clean and normalize gene expression data.

    Steps (aligned with WEKA):
    - Drop rows with missing values
    - Min-max normalize each gene column to [0, 1]
    """
    working = df.copy()
    label_col = label_col or _detect_label_column(working)
    gene_cols = _gene_columns(working)

    if not gene_cols:
        raise ValueError("No Gene_* columns found in the uploaded file.")

    missing_before = int(working.isna().sum().sum())
    rows_before = len(working)

    working = working.dropna(how="any")
    rows_dropped = rows_before - len(working)

    numeric_genes = working[gene_cols].apply(pd.to_numeric, errors="coerce")
    bad_numeric = int(numeric_genes.isna().sum().sum())
    if bad_numeric:
        working[gene_cols] = numeric_genes
        working = working.dropna(subset=gene_cols)
        rows_dropped = rows_before - len(working)

    cleaned_genes = working[gene_cols].copy()
    col_min = cleaned_genes.min()
    col_max = cleaned_genes.max()
    span = col_max - col_min
    span = span.replace(0, np.nan)
    cleaned_genes = (cleaned_genes - col_min) / span
    cleaned_genes = cleaned_genes.fillna(0.0)

    cleaned = pd.concat([cleaned_genes, working[[label_col]].reset_index(drop=True)], axis=1)

    summary = {
        "samples": len(cleaned),
        "genes": len(gene_cols),
        "label_column": label_col,
        "diagnosis_counts": cleaned[label_col].value_counts().to_dict(),
        "missing_values_removed": missing_before + bad_numeric,
        "rows_dropped": rows_dropped,
        "normalization": "Min-max per gene column (WEKA Normalize: scale=1, translate=0)",
    }
    return cleaned, summary


def top20_genes(cleaned: pd.DataFrame, label_col: str | None = None, n_genes: int = 20) -> tuple[pd.DataFrame, list[str]]:
    """
    Select top genes by ANOVA F-score and return long-format expression table.

    Matches WEKA AttributeEval + Ranker behavior used for leukemia_top20_genes.csv.
    """
    label_col = label_col or _detect_label_column(cleaned)
    gene_cols = _gene_columns(cleaned)

    f_scores, _ = f_classif(cleaned[gene_cols], cleaned[label_col])
    ranking = pd.Series(f_scores, index=gene_cols).sort_values(ascending=False)
    selected = list(ranking.head(n_genes).index)

    melted_rows: list[dict] = []
    for sample_idx, row in cleaned.iterrows():
        sample_id = sample_idx + 1
        diagnosis = row[label_col]
        for gene in selected:
            melted_rows.append(
                {
                    "Sample_ID": sample_id,
                    "Diagnosis": diagnosis,
                    "Gene": gene,
                    "Expression": float(row[gene]),
                }
            )

    return pd.DataFrame(melted_rows), selected


def to_csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")


def run_kmeans(cleaned_df: pd.DataFrame, n_clusters: int = 3, seed: int = 42) -> tuple[pd.DataFrame, dict]:
    """
    Run K-Means clustering using scikit-learn and return cluster assignments.
    """
    from sklearn.cluster import KMeans
    gene_cols = _gene_columns(cleaned_df)
    
    kmeans = KMeans(n_clusters=n_clusters, random_state=seed, n_init=10)
    clusters = kmeans.fit_predict(cleaned_df[gene_cols])
    
    cluster_labels = [f"cluster{c}" for c in clusters]
    label_col = _detect_label_column(cleaned_df)
    
    result = pd.DataFrame({
        "Instance_number": cleaned_df.index,
        "Diagnosis": cleaned_df[label_col],
        "Cluster": cluster_labels
    })
    
    counts = pd.Series(cluster_labels).value_counts().to_dict()
    return result, counts

