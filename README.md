# 🧬 ipyrad → IQ‑TREE: Ambystoma RADseq (single‑end, reference‑based)

This repository contains a complete workflow to:
1) run **ipyrad** (pre‑QC + steps 1–7),
2) **check output quality**, and
3) use the **PHYLIP** alignment to infer a **maximum‑likelihood tree in IQ‑TREE**.

Designed for *Ambystoma barbouri* & *A. texanum* single‑end RADseq, reference‑based assembly.

---

## 📦 Environment
```bash
mamba env create -f env/environment.yml
conda activate ipyrad-iqtree
