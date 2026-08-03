@echo off
setlocal
cd /d "%~dp0"

echo.
echo === Validate model against cymbal WAV recordings ===
echo Working directory: %CD%
echo.

REM Usage:
REM   run_validate.bat
REM   run_validate.bat "C:\path\to\wav\folder"
REM   run_validate.bat "C:\path\to\wav\folder" cymbal_46cm_medium
REM   run_validate.bat "C:\path\to\wav\folder" cymbal_46cm_medium "C:\path\to\report_dir"

set "WAV_DIR=%~1"
set "INSTRUMENT=%~2"
set "OUT_DIR=%~3"

if "%WAV_DIR%"=="" set "WAV_DIR=%~dp0wavs"
if "%INSTRUMENT%"=="" set "INSTRUMENT=cymbal_46cm_medium"
if "%OUT_DIR%"=="" set "OUT_DIR=%~dp0validation_out"

if not exist "%WAV_DIR%" (
  echo [FAIL] WAV folder not found: "%WAV_DIR%"
  echo.
  echo Put single-stroke cymbal WAVs in that folder, or pass a path:
  echo   run_validate.bat "D:\samples\cymbal"
  echo.
  pause
  exit /b 1
)

echo WAV dir   : %WAV_DIR%
echo Instrument: %INSTRUMENT%
echo Report dir: %OUT_DIR%
echo.

python validate_against_recordings.py --wav-dir "%WAV_DIR%" --instrument %INSTRUMENT% --out "%OUT_DIR%"
set ERR=%ERRORLEVEL%

echo.
if %ERR% neq 0 (
  echo [FAIL] validate_against_recordings.py exited with code %ERR%
) else (
  echo [OK] See "%OUT_DIR%\validation_report.md"
)
echo.
pause
exit /b %ERR%
