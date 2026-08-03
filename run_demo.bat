@echo off
setlocal
cd /d "%~dp0"

REM Compatibility wrapper: former run_demo.bat now launches NonTunPerc.
REM   run_demo.bat           -> GUI
REM   run_demo.bat --cli     -> headless pipeline (old demo behaviour)

call "%~dp0run_nontunperc.bat" %*
exit /b %ERRORLEVEL%
