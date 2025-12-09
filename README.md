Below is the **clean, final README.md** you can **directly paste into GitHub**.
I have extracted only the **important files from your folder** and ignored `.git/`, cache, logs, etc.

---

#  AI-Driven eDNA Classification & Biodiversity Analysis

### **Deep-Sea eDNA Taxonomy Identification | Streamlit Web App**

 **Live App:** [https://e-dna-classification.streamlit.app/](https://e-dna-classification.streamlit.app/)

This project provides an AI-powered, alignment-free pipeline to classify environmental DNA (eDNA), detect novel organisms, and compute biodiversity metrics using deep learning, FAISS vector search, and clustering methods.

---

##  Project Structure (Important Files Only)

```
E-DNA/
│
├── app.py                  # Main Streamlit app
├── app_v2.py               # Updated/experimental UI version
│
├── ssu_faiss.index         # FAISS index for fast vector search
├── ssu_vectorizer.pkl      # Vectorizer for embedding sequences
├── ssu_ids.pkl             # Taxonomy ID mapping
├── SSU_taxonomy.txt        # Taxonomy reference file
│
├── README.md               # Project documentation
└── .gitattributes
```

---

##  Features

###  **AI-Based Taxonomy Classification**

* Uses sequence embeddings
* Alignment-free (much faster than BLAST/QIIME)

###  **Fast FAISS Vector Search**

* Efficient large-scale taxonomy lookup
* GPU accelerated when available

###  **Novel Species Detection**

* HDBSCAN clustering highlights unknown/rare taxa

###  **Biodiversity Metrics**

Automatically computes:

* Shannon Index
* Simpson Index
* Chao1 species richness
* Abundance estimates

###  **Interactive Web App**

* Upload FASTA/CSV sequences
* Real-time predictions
* Visual biodiversity reports

---

##  Tech Stack

**Machine Learning:** PyTorch · Scikit-learn · NumPy · Pandas
**Vector Search:** FAISS
**Clustering:** HDBSCAN
**Web App:** Streamlit
**Data:** SSU rRNA taxonomy dataset

---

##  Getting Started

### 1️. Clone the Repository

```bash
git clone https://github.com/your-username/E-DNA.git
cd E-DNA
```

### 2️. Install Dependencies

```bash
pip install -r requirements.txt
```

### 3️. Run the App

```bash
streamlit run app.py
```

---

##  Input Format

Supports:

* `.fasta` sequences
* `.txt` raw reads
* `.csv` with sequence columns

---

##  Outputs

* Taxonomic prediction for each sequence
* Classification confidence
* Novelty alerts (unknown clusters)
* Biodiversity plots:

  * Taxa distribution
  * Cluster visualization
  * Diversity indices

---

##  Example Workflow

1. Upload your eDNA FASTA/CSV file
2. System cleans & embeds sequences
3. Performs vector search using FAISS
4. Groups unknowns using HDBSCAN
5. Generates biodiversity insights + visualizations

---

##  References

* DNABERT (Ji et al., 2021)
* SILVA SSU rRNA database
* HDBSCAN clustering
* FAISS similarity search

---

## 👥 Team

**Team: CixCodeCrushers**
Smart India Hackathon 2025
Problem Statement ID: *SIH25042*


