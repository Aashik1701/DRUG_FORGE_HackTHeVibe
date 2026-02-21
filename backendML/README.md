# backendML (Legacy Research Workspace)

This folder contains the older Flask-era research assets and notebooks used during model experimentation.

## Status

- `backendML` is **legacy/reference**.
- The active production API is in [../backend](../backend).
- Frontend integration should target the FastAPI service, not this folder.

## What is inside

```text
backendML/
├── app.py                         # Legacy Flask API (reference)
├── requirements.txt               # Legacy Python deps
├── ADMET Properties/              # Datasets, notebooks, saved models
├── Drug Target Binding Score/     # Binding model research assets
├── Molecular Docking/             # Docking notebooks and molecule files
└── Target Identification/         # Target-identification experiments
```

## When to use this folder

Use `backendML` only for:

- Reviewing historical experiments
- Retraining/rebuilding models
- Inspecting source datasets and notebooks

Do **not** use it as the primary API during hackathon demos.

## Active backend

Use [../backend](../backend) for:

- Running API locally (`uvicorn`)
- Swagger docs (`/docs`)
- Prediction endpoints (`/predict/*`)
- Health/model checks (`/health`, `/models`)

See root [README.md](../README.md) for full project flow.
