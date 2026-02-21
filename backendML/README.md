# 🧪 backendML — Model Research & Training

This folder contains **research assets, training notebooks, and model development code** for the AI models powering DrugForge.

---

## 🎯 Quick Answer: Which Backend Should I Use?

| Use Case | Folder | Command |
|----------|--------|---------|
| **🚀 Running the live API** | `../backend/` | `uvicorn main:app --reload` |
| **🔬 Training/retraining models** | `backendML/` | `python train_all_models.py` |
| **📚 Exploring model notebooks** | `backendML/` | Jupyter notebooks (see below) |
| **👨‍💻 Contributing new models** | `backendML/` | Add to `train_all_models.py` |

---

## 📁 Folder Structure

```
backendML/
│
├── 📄 README.md                    ← You are here
├── 📄 requirements.txt             ← Model training dependencies
├── 🐍 app.py                       ← Legacy Flask API (reference only)
├── 🔨 train_all_models.py          ← MAIN: Train all 9 models at once
│
├── 📊 ADMET Properties/            ← Absorption, Distribution, Metabolism, Excretion, Toxicity
│   ├── Datasets/
│   │   ├── bbb_martins.tab         ← Blood-Brain Barrier penetration data
│   │   ├── cyp3a4_veith.csv        ← CYP3A4 inhibition data
│   │   ├── half_life_obach.csv     ← Drug elimination half-life data
│   │   ├── herg_karim.tab          ← Heart toxicity (hERG) data
│   │   └── solubility_aqsoldb.tab  ← Water solubility data
│   └── Notebook/
│       ├── BBBP.ipynb              ← Train BBB model
│       ├── CYP P450 3A4 Inhibition.ipynb ← Train CYP3A4 model
│       ├── Solubility.ipynb        ← Train solubility model
│       ├── Toxicity.ipynb          ← Train toxicity models (HEPG2, BBBP-Tox)
│       └── Excretion.ipynb         ← Train half-life model
│
├── 🎯 Drug Target Binding Score/   ← Protein-ligand binding affinity
│   ├── Binding_Score.ipynb         ← Train binding prediction model
│   └── Datasets/                   ← (Binding training data)
│
├── 🧬 Molecular Docking/           ← 3D structure generation & docking
│   ├── molecular_docking.ipynb     ← Reference: AutoDock workflow
│   ├── sample.smi                  ← Example SMILES strings
│   └── *.mol2                      ← Example docked structures
│
└── 🎪 Target Identification/       ← Drug target prediction
    ├── Datasets/
    │   ├── ACE-2.csv               ← ACE2 receptor binding data
    │   └── COX-2.csv               ← COX-2 enzyme binding data
    └── Notebook/
        ├── ace2.ipynb              ← Train ACE2 binding model
        ├── cox2.ipynb              ← Train COX-2 binding model
        └── hepg2.ipynb             ← Train HepG2 toxicity model
```

---

## 🚀 How to Train Models

### Option 1: Train All Models at Once (Recommended)

```bash
cd backendML
python3 -m venv venv
source venv/bin/activate                    # Windows: venv\Scripts\activate
pip install -r requirements.txt
python train_all_models.py
```

**Output:** Trained `.pkl` files saved to `../backend/models/`

**Time:** ~5-10 minutes (depends on dataset sizes)

### Option 2: Train Individual Models (Jupyter Notebooks)

```bash
cd backendML
jupyter notebook
```

Then open any notebook in `ADMET Properties/Notebook/` or `Target Identification/Notebook/`

**Why:** Experiment with hyperparameters, visualize training curves, debug data

---

## 📊 Datasets Reference

| Dataset | File | Rows | Task | Model |
|---------|------|------|------|-------|
| BBB Penetration | `bbb_martins.tab` | 2,000+ | Binary classification | BBBP|
| CYP3A4 Inhibition | `cyp3a4_veith.csv` | 13,000+ | Binary classification | CYP3A4 |
| Water Solubility | `solubility_aqsoldb.tab` | 9,000+ | Regression | Solubility |
| Drug Toxicity | `herg_karim.tab` | 5,000+ | Binary classification | HEPG2, Toxicity |
| Drug Half-Life | `half_life_obach.csv` | 600+ | Regression | Half-Life |
| ACE2 Binding | `ACE-2.csv` | 100+ | Classification | ACE2 |
| COX-2 Binding | `COX-2.csv` | 400+ | Classification | COX-2 |

**Data Sources:**
- ChEMBL (EMBL-EBI)
- PubChem (NIH)
- MoleculeNet (Stanford)
- Literature research papers

---

## 🔧 For Contributors: Adding New Models

### 1. Add Your Dataset
```bash
mkdir -p backendML/New_Property/Datasets
# Place your .csv or .tab file here
```

### 2. Create a Training Notebook
```bash
# Create: backendML/New_Property/Notebook/my_model.ipynb
# Use existing notebooks as templates (e.g., BBBP.ipynb)
```

### 3. Update `train_all_models.py`
Add your model to the main training script:
```python
# In backendML/train_all_models.py
import your_notebook_module

def train_new_property():
    # Your training code
    model = ...
    joblib.dump(model, '../backend/models/new_property_model.pkl')
```

### 4. Update Backend Router
Create `../backend/routers/new_property.py`:
```python
@router.post("/predict/new-property")
def predict_new_property(payload: MoleculePayload):
    # Prediction logic
    pass
```

### 5. Test Locally
```bash
cd ../backend
python -m uvicorn main:app --reload --port 5001
# Test: curl http://localhost:5001/predict/new-property
```

---

## 🧬 Model Types & Algorithms

| Model | Algorithm | Type | Accuracy |
|-------|-----------|------|----------|
| **BBB Penetration** | Random Forest | Classification | 85-90% |
| **Solubility** | Gradient Boosting | Regression | MAE ~0.9 |
| **CYP3A4 Inhibition** | SVM | Classification | 80-85% |
| **Toxicity (HepG2)** | Neural Network | Classification | 90%+ |
| **ACE2 Binding** | Random Forest | Classification | 88%+ |
| **COX-2 Binding** | Gradient Boosting | Classification | 85%+ |
| **Half-Life** | Random Forest | Regression | MAE ~0.5 hrs |

**Feature Engineering:**
- RDKit molecular descriptors (200+ features)
- Fingerprints (ECFP, Morgan)
- Topological properties
- Electrostatic properties

---

## ⚠️ Important Notes

### ✅ DO
- Use this folder for model research & development
- Retrain models with new datasets
- Experiment with hyperparameters in notebooks
- Cite data sources in your changes

### ❌ DON'T
- Treat this as the production API
- Forget to export trained models to `../backend/models/`
- Skip testing after retraining
- Commit large dataset files (>10MB) without `.gitignore`

---

## 🔄 Workflow

### For Researchers
```
1. Explore data in Jupyter notebook
2. Train model with different hyperparameters
3. Evaluate performance (accuracy, MAE, etc.)
4. Export as .pkl file → ../backend/models/
5. Test in FastAPI → http://localhost:5001/predict/
```

### For Hackathon Competitors
```
1. The models are already trained (see ../backend/)
2. Use them via /predict/* endpoints
3. (If you want to improve: retrain here, then deploy to backend)
```

---

## 🔗 Related Files

- **Active Backend API:** [../backend/README.md](../backend/README.md)
- **Project Overview:** [../README.md](../README.md)
- **FastAPI Routes:** [../backend/routers/](../backend/routers/)
- **Model Loader:** [../backend/utils/model_loader.py](../backend/utils/model_loader.py)

---

## 📞 Questions?

- **"How do I retrain a model?"** → Run `train_all_models.py` after updating datasets
- **"Can I add a new prediction model?"** → Yes! Follow the "For Contributors" section above
- **"Where are the pre-trained models?"** → [../backend/models/](../backend/models/)
- **"Does the API use updated models automatically?"** → Yes, restart the API server

---

<div align="center">

**backendML** = The research & training lab  
**backend** = The production API lab  
Use the right tool for the right job! 🧪🚀

</div>
