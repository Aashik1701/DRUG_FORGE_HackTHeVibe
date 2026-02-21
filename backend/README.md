# DrugForge Backend (FastAPI)

Active backend service for DrugForge predictions.

## Quick Start

```bash
cd backend
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python -m uvicorn main:app --reload
```

- API root: http://localhost:8000/
- Swagger docs: http://localhost:8000/docs
- Health: http://localhost:8000/health
- Models: http://localhost:8000/models

## Main Endpoints

### Prediction endpoints (`POST /predict/*`)

- `/predict/solubility`
- `/predict/bbbp`
- `/predict/cyp3a4`
- `/predict/toxicity`
- `/predict/binding-score`
- `/predict/cox2`
- `/predict/hepg2`
- `/predict/ace2`
- `/predict/half-life`

Example request body:

```json
{
  "smiles": "CC(C)Cc1ccc(cc1)C(C)C"
}
```

### Utility endpoint

- `POST /utils/generate-3d` (returns 3D MOL block from SMILES)

## Environment Variables

- `MODELS_PATH` (optional): custom model directory
- `FRONTEND_URL` (optional): additional CORS origin
- `SUPABASE_URL` (optional): enable prediction persistence
- `SUPABASE_KEY` (optional): Supabase API key

If Supabase variables are missing, predictions still work; only persistence is skipped.

## Folder Overview

```text
backend/
├── main.py               # FastAPI app entry point
├── routers/              # Prediction and utility routers
├── schemas/              # Pydantic request/response models
├── services/             # DB service (Supabase)
├── utils/                # Model loader + RDKit helpers
├── models/               # Serialized ML models (.pkl)
└── tests/                # Smoke tests
```

## Testing

```bash
cd backend
source venv/bin/activate
python -m pytest tests/ -v
```
