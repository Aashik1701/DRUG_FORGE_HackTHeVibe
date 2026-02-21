# backendML Notes (Legacy)

This file documents the role of `backendML` after migration to FastAPI.

## Important

- `backendML` is kept for experiments and archived workflows.
- The currently maintained backend API is in `backend/`.

## Legacy run (only if you explicitly need it)

```bash
cd backendML
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python app.py
```

## Recommended run (active stack)

```bash
cd backend
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
python -m uvicorn main:app --reload
```

## Why both folders exist

- `backendML/`: historical model-building workspace
- `backend/`: production-ready API layer used by frontend
