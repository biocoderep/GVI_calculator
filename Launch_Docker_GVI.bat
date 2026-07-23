@echo off
echo ========================================================
echo GVI Calculator - Docker Edition Launcher
echo ========================================================
echo.
echo Ensuring Docker backend is running...
if not exist jobs.db type nul > jobs.db

set "GPU_FOUND=0"
nvidia-smi >nul 2>&1
if %ERRORLEVEL% equ 0 set "GPU_FOUND=1"

if exist "C:\Program Files\NVIDIA Corporation\NVSMI\nvidia-smi.exe" set "GPU_FOUND=1"
if exist "C:\Windows\System32\nvidia-smi.exe" set "GPU_FOUND=1"

if "%GPU_FOUND%"=="1" (
    echo [HARDWARE DETECTION] NVIDIA GPU Detected! Accelerating with GPU capabilities...
    docker-compose -f docker-compose.yml -f docker-compose.gpu.yml up --build -d
    if %ERRORLEVEL% neq 0 (
        echo [ERROR] Docker failed to mount the GPU ^(WSL driver issue^). Falling back to standard CPU mode...
        docker-compose up -d
    )
) else (
    echo [HARDWARE DETECTION] No NVIDIA GPU found. Launching standard CPU mode...
    docker-compose up --build -d
)

echo.
echo Launching your Web Dashboard!
start http://localhost:5000

echo.
echo The backend server is running safely in Docker.
echo You can close this black window at any time, the server will keep running in Docker!
echo.
pause
