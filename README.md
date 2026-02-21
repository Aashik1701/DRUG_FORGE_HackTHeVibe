# DrugForge: AI-Powered Drug Discovery Platform


---

## 🎯 What is DrugForge?

DrugForge rapidly predicts molecular properties and target interactions using machine learning. Scientists upload SMILES strings and get instant predictions on solubility, toxicity, half-life, BBB permeability, and more.

**Key Features:**
- 🔬 **Lab Bench** — Enter SMILES, get 9 simultaneous ML predictions
- 📊 **Batch Processor** — Upload CSV, predict 100s of molecules, export results
- 📱 **Professional UI** — Glass-morphism design, dark mode, responsive
- ⚡ **FastAPI Backend** — Async predictions with lazy-loaded models
- 🚀 **Production Ready** — Deploy to Vercel (frontend) + Render (backend)

---

## 🚀 Quick Start

### Prerequisites
- Python 3.10+
- Node.js 18+

### Backend Setup
```bash
cd backend
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python -m uvicorn main:app --reload
# API docs: http://localhost:8000/docs
```

### Frontend Setup (New Terminal)
```bash
npm install
npm run dev
# Frontend: http://localhost:5173
```

### Test
1. Open http://localhost:5173
2. Sign in (any email/password, 6+ chars)
3. Go to **Lab Bench**
4. Enter SMILES: `CC(C)Cc1ccc(cc1)C(C)C`
5. Click **Run All Models** → See predictions

---

## 🏗️ Architecture

```
Frontend (React + Vite)
  ↓ API calls (Axios)
FastAPI Backend (9 ML endpoints)
  ↓ Feature extraction (RDKit)
ML Models (scikit-learn, lazy-loaded)
  ↓ Optional: Save to DB (Supabase)
```

---

## 📋 Tech Stack

**Frontend:**
- React 18, Vite, Tailwind CSS, Framer Motion, RDKit.js, 3Dmol.js

**Backend:**
- FastAPI, Uvicorn, scikit-learn, RDKit, Pydantic


---


## 🎯 Main Features

### Lab Bench (Single Molecule Analysis)
- Enter SMILES string
- Run 9 simultaneous models
- See predictions with interpretations
- 2D/3D molecular visualization

### Batch Processor
- Upload CSV (auto-detects SMILES column)
- Select models to run
- Progress tracking
- Export results to CSV

### Prediction Models
- Solubility (logS)
- Toxicity (binary)
- BBBP (blood-brain barrier)
- CYP3A4 inhibition
- Half-life
- Binding score
- COX2, HepG2, ACE2

---

## 🔐 Security

- No user data stored (mock auth)
- HTTPS ready (Vercel + Render)
- CORS configured
- Secrets in .env.example template
