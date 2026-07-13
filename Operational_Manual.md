# GVI Calculator - Operational Manual

This manual provides detailed instructions on how to install, configure, and operate the Genomic Vulnerability Index (GVI) Calculator. The pipeline relies on heavy bioinformatics binaries (IQ-TREE 2, BEAST 2, PhiPack) to calculate advanced phylogenetic parameters.

You can run the application in two modes:
1. **Docker Mode (Recommended)**: Best for users who want a one-click setup without manual installations.
2. **Native / Non-Docker Mode**: Best for users who want to run the calculations natively on their hardware or debug the underlying scripts. For Windows users, this mode utilizes the Windows Subsystem for Linux (WSL) to execute Linux binaries natively.

---

## 🐳 Option 1: Docker Mode (Recommended)

In this mode, all dependencies, binaries, and environments are pre-configured inside an isolated container.

### Prerequisites
- **Docker Desktop** installed and running on your system.
- If you are on Windows, ensure Docker Desktop is configured to use the **WSL 2 backend** (this is usually the default).

### Installation & Execution
1. Open the project folder (`GVI_index_calculator`).
2. Double-click the **`Launch_Docker_GVI.bat`** script.
3. **What happens under the hood:**
   - The script runs `docker-compose up -d`.
   - Docker will build the image using the provided `Dockerfile`. It automatically downloads and configures a minimal Linux environment, installs Python, IQ-TREE 2, BEAST 2, and PhiPack directly into the container.
   - It mounts your `uploads/` folder and `jobs.db` database so your data persists even if the container stops.
   - The Flask web server starts on port `5000`.
4. The script will automatically open your default web browser to `http://localhost:5000`.
5. **To stop the server:** Open Docker Desktop and stop the container, or run `docker-compose down` in your terminal.

---

## 💻 Option 2: Native / Non-Docker Mode

In this mode, the Flask web server and the bioinformatics pipeline run directly on your host operating system. For Windows users, the web server runs in Windows, but it uses **WSL (Windows Subsystem for Linux)** to execute the heavy phylogenetic binaries (IQ-TREE, BEAST, Phi).

### Prerequisites
- **Python 3.10+** installed on Windows.
- **Windows Subsystem for Linux (WSL)** enabled and an Ubuntu/Debian distribution installed.
- **Java** installed inside WSL (required for BEAST 2).

### Step-by-Step Installation

#### 1. Install Python Dependencies
Open a standard Windows Command Prompt or PowerShell in the project directory and run:
```cmd
pip install -r requirements.txt
```

#### 2. Install IQ-TREE 2 and PhiPack (inside WSL)
Open your **WSL Terminal** (Ubuntu/Debian) and navigate to the project directory (e.g., `/mnt/e/Inlfuenza/avian_influenza/GVI_index_calculator/`). Run the provided shell script to download and install IQ-TREE and PhiPack:
```bash
chmod +x install_heavy_tools.sh
./install_heavy_tools.sh
```
*Note: The script installs IQ-TREE as `iqtree`. To ensure the Python backend finds it, create a symbolic link named `iqtree2`:*
```bash
sudo ln -s /usr/local/bin/iqtree /usr/local/bin/iqtree2
```

#### 3. Install BEAST 2 (inside WSL)
BEAST 2 requires Java to run. Open your **WSL Terminal** and execute the following:
```bash
# Install Java
sudo apt-get update
sudo apt-get install default-jre -y

# Navigate to the tools directory and extract BEAST
cd /mnt/e/Inlfuenza/avian_influenza/GVI_index_calculator/tools
tar -xzf BEAST.v2.7.6.Linux.x86.tgz

# Link the beast executable to your system path
sudo ln -s /mnt/e/Inlfuenza/avian_influenza/GVI_index_calculator/tools/beast/bin/beast /usr/local/bin/beast
```

### Execution
1. Open the project folder in Windows Explorer.
2. Double-click the **`Launch_GVI_Web_App.bat`** script.
3. **What happens under the hood:**
   - The script launches a WSL bash shell, navigates to the project directory, installs any missing python packages, and runs `python3 app.py`.
   - The Flask server boots up and binds to port `5000`.
   - When you upload a sequence and start a job, the backend `bio_engine.py` will dynamically pass commands to WSL (e.g., `wsl iqtree2 -s ...`) to perform the heavy number crunching.
4. The dashboard will automatically open at `http://localhost:5000`.
5. **To stop the server:** Close the black command prompt window that opened when you launched the `.bat` file.

---

## ⚙️ How the Pipeline Works

Regardless of which mode you run, the internal pipeline executes the following sequence:

1. **Alignment & Preprocessing**: Fasta sequences are chunked and aligned if necessary.
2. **Genetic Diversity**: Calculates Tajima's Pi (π).
3. **Codon Adaptation Index (CAI)**: Calculates species-specific adaptation (Sharp & Li).
4. **Selection Pressure (dN/dS)**: Estimates synonymous vs. non-synonymous mutations.
5. **ModelFinder**: IQ-TREE determines the optimal substitution model (e.g., GTR+I+G).
6. **BEAST 2 (Primary Clock & Recombination)**: Runs MCMC chains to calculate rigorous Recombination rates and the primary Evolutionary Rate.
7. **IQ-TREE LSD2 (Secondary Clock)**: Constructs a Maximum-Likelihood tree and calculates a secondary Evolutionary Rate (used as a fallback if BEAST fails to converge).
8. **BDSKY (Reproductive Number)**: Runs a Birth-Death Skyline model in BEAST 2 to estimate the Effective Reproductive Number (Re).
9. **Final GVI Calculation**: Normalizes the 9 parameters into a single `0.0` - `1.0` index score.

### Troubleshooting (Native Mode)
- **"Command not found: beast"**: Ensure you ran the `ln -s` command for BEAST in your WSL terminal, and that the path strictly points to the `bin/beast` executable.
- **"IQ-TREE execution failed"**: Ensure you created the `iqtree2` symlink. The python backend strictly looks for `iqtree2`.
- **BEAST crashes immediately**: Usually caused by missing Java. Ensure `default-jre` is installed in WSL.
