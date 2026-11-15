@echo off
REM Double-click this file to start NovaScope with real SPHEREx data.
REM Keep the black window open while you use NovaScope. Close it (or press Ctrl+C) to stop.
title NovaScope
cd /d "%~dp0"
if not exist "backend\app.py" (
  echo Put this file in the novascope folder, next to the "backend" and "frontend" folders.
  pause
  exit /b 1
)
where py >nul 2>nul
if %errorlevel%==0 (
  py backend\app.py
) else (
  python backend\app.py
)
echo.
echo NovaScope has stopped.
pause
