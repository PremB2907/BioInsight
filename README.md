# 🧬 Leukemia Bio-Intelligence SLA Dashboard

An interactive, high-fidelity web application built with **Flask + Vanilla HTML/CSS/JS** and **Chart.js** to preprocess, select, and cluster leukemia gene expression profiles. This portal mirrors the exact analytical pages and layout logic of the Power BI report (`DMBI SLA.pbix`) and aligns with the WEKA data mining algorithms.

---

## 🚀 Key Features

* **Advanced Normalization (WEKA SimpleKMeans standard)**: Cleans datasets of missing rows and executes Min-Max Normalization to scale expression levels strictly between `[0, 1]`.
* **ANOVA F-Score Attribute Selector**: Reduces dimensionality from 1,000 genes to the top 20 genes explaining the highest statistical variance between classes.
* **Interactive Permutation Matcher**: Evaluates unsupervised cluster geometries against clinical labels in real-time, displaying accuracy readouts for all $3! = 6$ possible mappings.
* **Biometric Boxplots & Gene Profiles**: Features a custom HTML5 canvas boxplot renderer displaying standard deviations and mean lines, alongside interactive line charts of patient scans.
* **Drag-and-Drop Preprocessing Pipeline**: Upload custom CSV/ARFF datasets on-the-fly, run the extraction pipeline, and export the processed datasets immediately.

---

## 📊 Scientific & Data Mining Workflow

```mermaid
graph TD
    A[Raw Gene Expression: 1000 Genes + Diagnosis] --> B[Min-Max Normalization to [0,1]]
    B --> C[ANOVA F-Score Feature Ranking]
    C --> D[Select Top 20 Genes]
    B --> E[Unsupervised SimpleKMeans Clustering K=3]
    E --> F[Interactive Permutation Alignments]
    F --> G[Classification Accuracies]
```

### 1. Preprocessing & Normalization
Eliminates scaling bias using **Min-Max scaling** on gene expression levels:
$$X_{norm} = \frac{X - X_{min}}{X_{max} - X_{min}}$$
This prevents extremely active gene probes from skewing distance evaluations.

### 2. Feature Selection (ANOVA F-score)
Computes the variance ratio to rank attributes by classification significance:
$$F = \frac{\text{Variance between diagnosis classes}}{\text{Variance within diagnosis classes}}$$
The top 20 features are selected to solve the *curse of dimensionality*.

### 3. Unsupervised Clustering
Partitions the 1000-gene normalized space into $K=3$ clusters. An interactive permuter resolves the arbitrary cluster-to-diagnosis associations, locating the optimal mapping of `cluster0 ➔ AML`, `cluster1 ➔ Healthy`, `cluster2 ➔ ALL` yielding **75.20% accuracy**.

---

## 🛠️ Installation & Setup

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/PremB2907/BioInsight.git
   cd BioInsight
   ```

2. **Install Dependencies**:
   Ensure you have Python 3.10+ installed. Install dependencies using:
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the Server**:
   Launch the Flask server:
   ```bash
   python app.py
   ```

4. **Access the Portal**:
   Open your browser and navigate to:
   👉 **http://localhost:5000**

---

## 📁 Repository Structure

```
├── app.py                     # Flask server exposing REST JSON APIs
├── processing.py              # Data cleaning, ANOVA selection, and KMeans logic
├── requirements.txt           # Python library dependencies
├── templates/
│   └── index.html             # Cyper-aesthetic single page UI, styles, and Chart.js code
├── DMBI SLA.pbix              # Underlying Power BI dashboard structure
├── cluster_assignments.csv     # Precalculated WEKA SimpleKMeans clusters
├── leukemia_top20_genes.csv   # Long-format ranked gene expression CSV
└── leukemia_gene_expression_cleaned.csv  # Min-Max normalized full cohort
```
