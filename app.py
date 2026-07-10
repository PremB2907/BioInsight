from flask import Flask, jsonify, request, render_template
import pandas as pd
import numpy as np
import io
import os

# Import our custom processing pipeline functions
from processing import clean_dataset, load_dataset, top20_genes, run_kmeans

app = Flask(__name__)

# Cache variables for local data files
CLEANED_PATH = "Data/leukemia_gene_expression_cleaned.csv"
TOP20_PATH = "Data/leukemia_top20_genes.csv"
CLUSTERS_PATH = "Data/cluster_assignments.csv"

# Pre-load data if available
def load_data():
    data = {}
    if os.path.exists(CLEANED_PATH):
        data["cleaned"] = pd.read_csv(CLEANED_PATH)
    else:
        data["cleaned"] = None

    if os.path.exists(TOP20_PATH):
        data["top20"] = pd.read_csv(TOP20_PATH)
    else:
        data["top20"] = None

    if os.path.exists(CLUSTERS_PATH):
        data["clusters"] = pd.read_csv(CLUSTERS_PATH)
    else:
        data["clusters"] = None
    return data

DATA_CACHE = load_data()

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/api/stats", methods=["GET"])
def get_stats():
    global DATA_CACHE
    cleaned_df = DATA_CACHE["cleaned"]
    clusters_df = DATA_CACHE["clusters"]
    top20_df = DATA_CACHE["top20"]
    
    if cleaned_df is None or clusters_df is None or top20_df is None:
        # Reload cache in case files were generated
        DATA_CACHE = load_data()
        cleaned_df = DATA_CACHE["cleaned"]
        clusters_df = DATA_CACHE["clusters"]
        top20_df = DATA_CACHE["top20"]

    if cleaned_df is None or clusters_df is None or top20_df is None:
        return jsonify({"error": "Local datasets not found. Please upload a dataset in the pipeline tab."}), 404

    total_samples = len(cleaned_df)
    gene_cols = [c for c in cleaned_df.columns if c.startswith("Gene_")]
    total_genes = len(gene_cols)
    
    diag_counts = cleaned_df["Diagnosis"].value_counts().to_dict()
    cluster_counts = clusters_df["Cluster"].value_counts().to_dict()
    
    # Generate Crosstab Matrix
    ct = pd.crosstab(clusters_df["Diagnosis"], clusters_df["Cluster"])
    crosstab_data = {index: row.to_dict() for index, row in ct.iterrows()}
    
    # Get top 20 genes list
    top20_genes_list = sorted(top20_df["Gene"].unique().tolist())
    
    # Send cluster assignments so JS can compute accuracy locally on-the-fly
    cluster_assignments = clusters_df.to_dict(orient="records")

    return jsonify({
        "total_samples": total_samples,
        "total_genes": total_genes,
        "diagnosis_counts": diag_counts,
        "cluster_counts": cluster_counts,
        "crosstab": crosstab_data,
        "top20_genes": top20_genes_list,
        "cluster_assignments": cluster_assignments
    })

@app.route("/api/gene-stats/<gene_name>", methods=["GET"])
def get_gene_stats(gene_name):
    cleaned_df = DATA_CACHE["cleaned"]
    top20_df = DATA_CACHE["top20"]
    
    if cleaned_df is None or top20_df is None:
        return jsonify({"error": "Datasets not loaded."}), 404

    if gene_name not in cleaned_df.columns:
        return jsonify({"error": f"Gene {gene_name} not found."}), 404

    # Calculate statistics grouped by Diagnosis
    stats = cleaned_df.groupby("Diagnosis")[gene_name].agg(["mean", "std", "min", "max", "count"]).to_dict(orient="index")
    
    # Extract expressions and diagnoses for the line graph
    expression_data = []
    for idx, row in cleaned_df.iterrows():
        expression_data.append({
            "sample_index": idx + 1,
            "diagnosis": row["Diagnosis"],
            "expression": float(row[gene_name])
        })
        
    return jsonify({
        "gene": gene_name,
        "stats": stats,
        "expressions": expression_data
    })

@app.route("/api/process", methods=["POST"])
def process_file():
    if "file" not in request.files:
        return jsonify({"error": "No file uploaded."}), 400
        
    uploaded_file = request.files["file"]
    if uploaded_file.filename == "":
        return jsonify({"error": "Empty filename."}), 400

    try:
        file_bytes = uploaded_file.read()
        raw_df = load_dataset(file_bytes, uploaded_file.filename)
        
        # Run clean pipeline
        cleaned_df, summary = clean_dataset(raw_df)
        
        # Run feature selection (ANOVA)
        top20_df, selected_genes = top20_genes(cleaned_df, summary["label_column"])
        
        # Run clustering (K-Means)
        cluster_assignments, cluster_counts = run_kmeans(cleaned_df)
        
        # Convert to serializable objects
        cleaned_data_sample = cleaned_df.head(10).to_dict(orient="records")
        top20_data = top20_df.to_dict(orient="records")
        cluster_data = cluster_assignments.to_dict(orient="records")
        
        # Cleaned dataset full export via CSV string
        csv_cleaned = cleaned_df.to_csv(index=False)
        csv_top20 = top20_df.to_csv(index=False)
        csv_clusters = cluster_assignments.to_csv(index=False)
        
        return jsonify({
            "success": True,
            "summary": summary,
            "selected_genes": selected_genes,
            "cluster_counts": cluster_counts,
            "cleaned_sample": cleaned_data_sample,
            "csv_cleaned": csv_cleaned,
            "csv_top20": csv_top20,
            "csv_clusters": csv_clusters
        })
        
    except Exception as e:
        return jsonify({"error": str(e)}), 500

if __name__ == "__main__":
    # Run the Flask app on port 5000
    app.run(host="0.0.0.0", port=5000, debug=True)
