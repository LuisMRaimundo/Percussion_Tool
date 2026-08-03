@echo off
setlocal
cd /d "%~dp0"

echo.
echo === NonTunPerc — Validate against recordings ===
echo Working directory: %CD%
echo.

REM Default: open GUI (browse sample folder; searches subfolders by default).
REM Accepts .wav / .aif / .aiff / .flac.
REM CLI:
REM   run_validate.bat --cli "C:\path\to\sample\folder"
REM   run_validate.bat --cli "C:\path\to\sample\folder" cymbal_46cm_medium
REM   run_validate.bat --cli "C:\path\to\sample\folder" cymbal_46cm_medium "C:\path\to\report_dir"

if "%~1"=="" (
  echo Mode: GUI
  python validate_against_recordings.py --gui
  exit /b %ERRORLEVEL%
)

if /I "%~1"=="--gui" (
  python validate_against_recordings.py --gui
  exit /b %ERRORLEVEL%
)

if /I "%~1"=="--cli" (
  set "WAV_DIR=%~2"
  set "INSTRUMENT=%~3"
  set "OUT_DIR=%~4"
) else (
  REM Backward compatible: first arg is wav folder
  set "WAV_DIR=%~1"
  set "INSTRUMENT=%~2"
  set "OUT_DIR=%~3"
)

if "%WAV_DIR%"=="" (
  echo [FAIL] Sample folder required for CLI mode.
  echo   run_validate.bat --cli "D:\samples\cymbal"
  pause
  exit /b 1
)
if "%INSTRUMENT%"=="" set "INSTRUMENT=cymbal_46cm_medium"
if "%OUT_DIR%"=="" set "OUT_DIR=%~dp0validation_out"

if not exist "%WAV_DIR%" (
  echo [FAIL] Sample folder not found: "%WAV_DIR%"
  pause
  exit /b 1
)

echo Mode: CLI
echo Sample dir: %WAV_DIR%
echo Recursive : yes (default; use python --no-recursive to disable)
echo Instrument: %INSTRUMENT%
echo Report dir: %OUT_DIR%
echo.

python validate_against_recordings.py --cli --wav-dir "%WAV_DIR%" --instrument %INSTRUMENT% --out "%OUT_DIR%"
set ERR=%ERRORLEVEL%

echo.
if %ERR% neq 0 (
  echo [FAIL] validation exited with code %ERR%
) else (
  echo [OK] See "%OUT_DIR%\validation_report.md"
)
echo.
pause
exit /b %ERR%
