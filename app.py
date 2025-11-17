# app.py
import streamlit as st
import pickle
import faiss
import numpy as np
from Bio import SeqIO
import io
import re
import pandas as pd
import os

# ---------------- CONFIG ----------------
FAISS_INDEX_PATH   = "ssu_faiss.index"        # adjust paths if needed
VECTORIZER_PATH    = "ssu_vectorizer.pkl"
IDS_PATH           = "ssu_ids.pkl"
TAX_PATH           = "ssu_taxonomy.pkl"
DEFAULT_KMER       = 6
DEFAULT_TOPK       = 5
DEFAULT_THRESHOLD  = 1.2   # distance threshold that we have decided 
# ----------------------------------------------------------------

st.set_page_config(page_title="eDNA SSU Classifier", layout="wide")

st.title(" eDNA SSU Classifier Prototype")
st.write("Paste a DNA sequence or upload a FASTA file. The system will classify it against the *SSU reference DB* and mark it as *Known* or *Novel candidate*.")

#  Define tokenizer before loading vectorizer
def split_tokens(s):
    return s.split()

# ---------------- Load artifacts ----------------
@st.cache_resource(ttl=3600)
def load_artifacts():
    if not os.path.exists(FAISS_INDEX_PATH):
        st.error(f"FAISS index not found at {FAISS_INDEX_PATH}")
        return None
    index = faiss.read_index(FAISS_INDEX_PATH)
    with open(VECTORIZER_PATH, "rb") as f:
        vectorizer = pickle.load(f)
    with open(IDS_PATH, "rb") as f:
        ids = pickle.load(f)
    # taxonomy is optional
    tax_map = {}
    if os.path.exists(TAX_PATH):
        try:
            with open(TAX_PATH, "rb") as f:
                tax_map = pickle.load(f)
        except Exception:
            tax_map = {}
    return {"index": index, "vectorizer": vectorizer, "ids": ids, "tax_map": tax_map}

art = load_artifacts()
if art is None:
    st.stop()

index = art["index"]
vectorizer = art["vectorizer"]
ids = art["ids"]
tax_map = art["tax_map"]

# ---------------- Helpers ----------------
def clean_seq(s):
    return re.sub(r'[^ACGT]', '', s.upper())

def seq_to_kmers(seq, k=6):
    seq = seq.upper()
    if len(seq) < k:
        return ""
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    kmers = [km for km in kmers if all(ch in "ACGT" for ch in km)]
    return " ".join(kmers)

def interpret_distance(dist, thresholds=(1.0, 1.5, 2.0)):
    if dist <= thresholds[0]:
        return "Very close — likely same species"
    if dist <= thresholds[1]:
        return "Close — likely same genus"
    if dist <= thresholds[2]:
        return "Distant — possible novel lineage"
    return "Very distant — novel/unknown"

global label

def classify_query(seq, kmer, top_k, threshold):
    kmers = seq_to_kmers(seq, k=kmer)
    if not kmers:
        return {"error": "Sequence too short or invalid."}
    emb = vectorizer.transform([kmers]).toarray().astype("float32")
    D, I = index.search(emb, top_k)
    neighbors = []
    for rank, (dist, idx) in enumerate(zip(D[0], I[0]), start=1):
        if idx < 0 or idx >= len(ids):
            continue
        ref_id = ids[idx]
        tax = tax_map.get(ref_id, "Unknown")
        neighbors.append({"rank": rank, "ref_id": ref_id, "distance": float(dist), "taxonomy": tax})
    best_dist = neighbors[0]["distance"] if neighbors else None
    label = "Novel candidate" if (best_dist is None or best_dist > threshold) else "Known (close match)"
    return {"label": label, "best_distance": best_dist, "neighbors": neighbors}

# ---------------- Sidebar controls ----------------
with st.sidebar:
    st.header("⚙ Parameters")
    kmer = st.number_input("k-mer length (k)", value=DEFAULT_KMER, min_value=3, max_value=12, step=1)
    top_k = st.number_input("Top matches", value=DEFAULT_TOPK, min_value=1, max_value=20, step=1)
    threshold = st.slider("Novelty threshold (distance)", min_value=0.0, max_value=5.0, value=float(DEFAULT_THRESHOLD), step=0.01)
    st.caption("Distance ≤ threshold → Known; > threshold → Novel candidate")

# ---------------- Input ----------------
st.subheader("Input Sequence")
col1, col2 = st.columns([2,1])
with col1:
    seq_input = st.text_area("Paste a DNA sequence ", height=150)
    uploaded = st.file_uploader("Or upload a FASTA file", type=["fa", "fasta", "fas", "txt"])
with col2:
    if st.button("Use example sequence"):
        seq_input = "AGCCAGGCATGTCTAGTACAAACCCCAAAGGGGGAAAACCGCGAAAGGCTCATTAAATCAGTTATGTTTCCTTTGATTGGAACCTTTTTTTTGATAACGGTGGTAATTCTAGAGCAAATA"
        st.experimental_rerun()

if uploaded:
    text = uploaded.getvalue().decode('utf-8')
    fasta_io = io.StringIO(text)
    recs = list(SeqIO.parse(fasta_io, "fasta"))
    if len(recs) > 0:
        seq_input = str(recs[0].seq)
        st.success(f"Loaded FASTA record: {recs[0].id}")

if not seq_input or seq_input.strip() == "":
    st.info("Paste a sequence or upload a FASTA to classify.")
    st.stop()

# ---------------- Run classification ----------------
seq_clean = clean_seq(seq_input)
if len(seq_clean) < kmer:
    st.error(f"Sequence too short after cleaning ({len(seq_clean)} bases). Need ≥ {kmer}.")
    st.stop()

with st.spinner("Running classification..."):
    res = classify_query(seq_clean, kmer, top_k, threshold)

if "error" in res:
    st.error(res["error"])
    st.stop()

# ---------------- Results ----------------
st.subheader("Result")
colA, colB = st.columns([2,1])
with colA:
    st.markdown(f"*Label:* {res['label']}")
    st.markdown(f"*Best distance:* {res['best_distance']:.4f}")
    st.markdown(f"*Interpretation:* {interpret_distance(res['best_distance'])}")
with colB:
    st.markdown("*Summary*")
    st.write(f"Seq length: {len(seq_clean)} bp")
    st.write(f"k-mer size: {kmer}")
    st.write(f"Top k: {top_k}")

if res["neighbors"]:
    df = pd.DataFrame(res["neighbors"])
    df = df.rename(columns={"ref_id":"Reference ID", "taxonomy":"Taxonomy", "distance":"Distance", "rank":"Rank"})
    st.dataframe(df)

st.markdown("---")
st.caption("Prototype using TF-IDF k-mer embeddings + FAISS. For production → replace with DNABERT embeddings + benchmark against QIIME2/DADA2.")