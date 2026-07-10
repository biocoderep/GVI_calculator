@echo off
echo ========================================================
echo GVI Calculator - Ultimate Edition (Web App Launcher)
echo ========================================================
echo.
echo Starting the backend Supercomputer (Flask)...
start /b wsl bash -i -c "cd /mnt/e/Inlfuenza/avian_influenza/GVI_index_calculator && pip install -r requirements.txt && python3 app.py"

echo Waiting for the server to boot up...
timeout /t 5 /nobreak > nul

echo Launching your Web Dashboard!
start http://localhost:5000

echo.
echo The Dashboard is now open in your web browser.
echo Do NOT close this black window until you are finished calculating!
echo To stop the server, just close this window.
pause
