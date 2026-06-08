@echo off
setlocal

set "SCRIPT_DIR=%~dp0"
set "PROJECT_ROOT=%SCRIPT_DIR%..\.."
set "VENV_PY=%PROJECT_ROOT%\.venv\Scripts\python.exe"

if exist "%VENV_PY%" (
    set "PYTHON_CMD=%VENV_PY%"
) else (
    set "PYTHON_CMD=python"
)

cd /d "%SCRIPT_DIR%"
echo Launching AnDi diffusion demo with: %PYTHON_CMD%
"%PYTHON_CMD%" "%SCRIPT_DIR%run.py"
