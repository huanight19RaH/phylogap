@echo off
title PhyloGap Studio Web Dashboard
echo ========================================================
echo       PhyloGap Studio - Genome Assembly Gap Reconstruction
echo       B10K Avian Calibrated Multimodal Engine
echo ========================================================
echo.
cd /d "%~dp0"

echo [+] Starting backend server on http://127.0.0.1:8085 ...
start "" "http://127.0.0.1:8085"
python -m uvicorn phylogap_ui.server:app --host 127.0.0.1 --port 8085 --reload

pause
