# Leukemia Gene Expression Analytics Dashboard

An interactive Power BI dashboard analyzing gene expression data across Acute Myeloid Leukemia (AML), Acute Lymphoblastic Leukemia (ALL), and Healthy control samples, combined with unsupervised clustering (K-Means via WEKA) and association rule mining to uncover patterns across a 20-gene expression panel.

## Overview

This project explores whether gene expression profiles can meaningfully distinguish between AML, ALL, and healthy samples. It combines exploratory data analysis, feature selection, unsupervised clustering, and association rule mining, with results visualized in a 4-page interactive Power BI dashboard.

## Dataset

- **Source:** [Leukemia Gene Expression Dataset (Kaggle)](https://www.kaggle.com/datasets/ziya07/gene-expression-dataset)
- **Samples:** 1,000 patient samples
- **Features:** 1,000 gene expression attributes (reduced to top 20 via feature selection)
- **Classes:**
  - AML — 409 samples (40.9%)
  - ALL — 380 samples (38%)
  - Healthy — 211 samples (21.1%)

## Tools & Technologies

| Tool | Purpose |
|---|---|
| **WEKA** | Feature selection (Information Gain), K-Means clustering, Apriori association rule mining |
| **Python** | Data preprocessing, reshaping (wide-to-long), cluster export handling |
| **Power BI** | Data modeling, DAX measures, interactive dashboard |
| **Excel** | Initial data inspection |

## Repository Structure

```
├── data/
│   ├── leukemia_gene_expression.csv           # Raw dataset (1000 samples x 1000 genes)
│   ├── leukemia_gene_expression_cleaned.csv   # Normalized dataset
│   ├── leukemia_top20_genes.csv               # Top 20 genes (wide format, post feature selection)
│   ├── leukemia_top20_genes_long.csv          # Unpivoted long format for Power BI
│   └── cluster_assignments.csv                # Sample-to-cluster mapping (from WEKA K-Means)
├── weka/
│   └── cluster_visualization.arff             # WEKA clustering output
├── dashboard/
│   └── leukemia_dashboard.pbix                # Power BI dashboard file
└── README.md
```

## Methodology

### 1. Data Preprocessing
- Missing value analysis and validation
- Normalization of gene expression values (0–1 scale)
- Feature selection using **Information Gain Attribute Evaluator** in WEKA to reduce 1,000 genes down to the top 20 most informative genes for classifying AML / ALL / Healthy

### 2. Data Reshaping (Python)
- The top-20-genes dataset was reshaped from wide format (one column per gene) to long format (`Sample_ID, Diagnosis, Gene, Expression`) to support proper Power BI data modeling — enabling gene-level slicers, distribution charts, and heatmaps without needing 20 separate visuals.

### 3. Exploratory Data Analysis
- Diagnosis-wise sample distribution
- Mean vs. median expression comparison across AML / ALL / Healthy
- Identification of top differentially expressed genes (by expression range across diagnosis groups)

### 4. Data Mining (WEKA)
- **K-Means Clustering (k=3, Euclidean distance):** Clustering was run independent of diagnosis labels to check whether natural groupings in the gene expression data align with known disease classes.
  - Cluster 0: 457 samples (45.7%)
  - Cluster 1: 94 samples (9.4%)
  - Cluster 2: 449 samples (44.9%)
  - Clusters showed partial alignment with actual diagnosis labels but did not achieve complete separation, suggesting the 20-gene subset captures meaningful but incomplete biological variation.
- **Association Rule Mining (Apriori):** Discovered 203 frequent itemsets and 28 large 2-itemsets after discretization. No rules met the confidence threshold — indicating weak but observable gene co-occurrence patterns rather than strong deterministic rules.
- **Outlier Detection:** Visualized via expression distribution spread per diagnosis group.

### 5. Data Modeling (Power BI)
A star schema was used:
- **Fact table:** `leukemia_top20_genes` (long format — Sample_ID, Diagnosis, Gene, Expression)
- **Dim_Diagnosis:** Diagnosis, SortOrder, ColorHex
- **Dim_Gene:** unique gene list
- **cluster_assignments:** Sample_ID, Cluster (joined 1:many to fact table)

**Key DAX measures:**
```DAX
Avg Expression = AVERAGE(leukemia_top20_genes[Expression])
Median Expression = MEDIAN(leukemia_top20_genes[Expression])
Sample Count = DISTINCTCOUNT(leukemia_top20_genes[Sample_ID])
```

## Dashboard Pages

| Page | Contents |
|---|---|
| **Executive Overview** | KPI cards (Total Samples, AML/ALL/Healthy counts), Diagnosis distribution donut, Mean vs Median expression chart, Top differentially expressed genes |
| **Gene Analysis** | Expression trend by diagnosis, Gene × Diagnosis heatmap matrix, interactive gene slicer, feature selection methodology notes |
| **Cluster Analysis** | Sample-level scatter plot colored by cluster, cluster size distribution, K-Means methodology notes |
| **Summary** | Key findings across dataset, gene analysis, clustering, and association rule mining, with recommendations |

## Key Findings

- Average and median expression values were closely aligned across all three diagnosis groups, indicating fairly symmetric distributions within the 20-gene panel.
- Genes with the largest expression range across diagnosis groups (strongest candidates for differentiation): **Gene_244, Gene_521, Gene_366, Gene_646, Gene_149**.
- K-Means clustering (k=3) partially aligned with actual diagnosis labels but did not fully separate the three classes — suggesting this 20-gene subset alone provides moderate, not complete, natural separation.
- Apriori association rule mining found frequent gene co-occurrence patterns but no rules strong enough to meet standard confidence thresholds.

## Limitations

- Analysis is restricted to the top 20 genes selected via Information Gain; the full 1,000-gene feature space was not clustered or mined directly.
- No hard performance benchmarks (e.g., precision/recall against known biomarkers) were available for cross-validation within the scope of this dataset.
- Dataset is Kaggle-sourced and should be treated as exploratory rather than clinically validated.

## Recommendations / Future Work

1. Validate top differentiator genes against known leukemia biomarker literature.
2. Re-run clustering on the full 1,000-gene feature set for stronger class separation.
3. Cross-validate K-Means clusters against Diagnosis labels using a confusion-matrix comparison.
4. Explore supervised classification models (e.g., Random Forest, SVM) as a complementary approach to unsupervised clustering.
5. Extend association rule mining with lower confidence thresholds to surface weaker but potentially informative patterns.

## How to Use

1. Clone this repository.
2. Open `dashboard/leukemia_dashboard.pbix` in Power BI Desktop.
3. Data source paths may need to be updated under **Transform Data → Data Source Settings** if CSVs are moved.
4. Use the gene slicer and page navigation buttons to explore each analysis page.

## Author

Prem Baraskar
