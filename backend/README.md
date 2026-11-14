# ASTRAEA Backend (migration starting point)

This folder contains the latest backend module implementations extracted from the supplied `NASA.ipynb`. Only the leading Jupyter `%%writefile` directives were removed; implementation logic has not yet been refactored or fully validated.

## Setup

From this `backend/` directory:

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

## Start the API

```bash
uvicorn main:app --reload
```

Open `http://127.0.0.1:8000/docs` to inspect the current API.

## Current endpoints from the notebook

- `GET /api/targets`
- `GET /api/frames?ra=...&dec=...`
- `GET /api/spectrum?ra=...&dec=...`

## Important

This is a migration baseline, not a production-ready backend. Run it locally, record any import/runtime errors, and do not claim scientific validation until the processing and outputs have been checked. The current notebook implementation still needs work on object-name resolution, image validation, detector-compatible epoch comparisons, serialization of non-finite values, error handling, caching policy, and automated/scientific tests.
