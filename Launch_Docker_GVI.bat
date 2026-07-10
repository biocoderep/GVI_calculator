@echo off
echo ========================================================
echo GVI Calculator - Docker Edition Launcher
echo ========================================================
echo.
echo Ensuring Docker backend is running...
docker-compose up -d

echo.
echo Launching your Web Dashboard!
start http://localhost:5000

echo.
echo The backend server is running safely in Docker.
echo You can close this black window at any time, the server will keep running in Docker!
echo.
pause
