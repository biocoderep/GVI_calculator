# GVI Calculator

![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![Flask](https://img.shields.io/badge/flask-%23000.svg?style=for-the-badge&logo=flask&logoColor=white)

**GVI Calculator** is a bioinformatics pipeline and interactive Web Dashboard designed to quantify the evolutionary risk of viral populations (originally designed for Avian Influenza). 

By automatically integrating advanced phylogenetic tools, the pipeline synthesizes 9 critical genomic parameters—including Evolutionary Rate, Effective Reproductive Number (Re), and Recombination—into a single normalized Genomic Vulnerability Index (GVI) on a 0.0 to 1.0 scale.

## Key Features
* **Full Phylogenetic Integration**: Automates execution of `IQ-TREE 2` and Headless `BEAST 2` natively.
* **Advanced Phylodynamics**: Automatically estimates molecular clock rates using `LSD2` and transmission dynamics using `BDSKY`.
* **Interactive UI**: A beautiful, real-time web dashboard to visualize genomic risk metrics over time.
* **100% Reproducible**: Fully Dockerized backend ensuring cross-platform stability without complex dependency management.

## Installation & Execution

### Option 1: Docker (Recommended)
This application is fully Dockerized for maximum cross-platform compatibility without needing to configure C++ or Java binaries manually.

1. Ensure you have [Docker Desktop installed](https://www.docker.com/products/docker-desktop/) and running.
2. Clone this repository and navigate into it:
   ```bash
   git clone https://github.com/biocoderep/GVI_calculator.git
   cd GVI_calculator
   ```
3. Run the following command:
   ```bash
   docker-compose up --build -d
   ```
4. Open your web browser and navigate to: [http://localhost:5000](http://localhost:5000)

### Option 2: Native Installation (No Docker)
If you cannot use Docker, you must install the dependencies manually on your operating system (Linux or macOS recommended; Windows users should use WSL).

1. **Install Python 3.10+**: Ensure Python and `pip` are installed.
2. **Clone the repository**:
   ```bash
   git clone https://github.com/biocoderep/GVI_calculator.git
   cd GVI_calculator
   ```
3. **Install IQ-TREE 2**: Download [IQ-TREE 2](http://www.iqtree.org/) and add the `iqtree2` binary to your system PATH.
4. **Install BEAST 2**: Download [BEAST 2](https://www.beast2.org/), ensure Java is installed, and add the `beast` binary to your system PATH.
5. **Install Python Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
5. **Run the Application**:
   ```bash
   python app.py
   ```
6. Open your web browser and navigate to: [http://localhost:5000](http://localhost:5000)

## Mathematical Parameters Calculated
The pipeline dynamically processes FASTA sequences and derives:
1. **Selection Pressure (dN/dS)**: Nei-Gojobori / Jukes-Cantor
2. **Codon Adaptation Index (CAI)**: Sharp & Li
3. **Evolutionary Rate**: ModelFinder + IQ-TREE 2 + LSD2
4. **Effective Reproductive Number (Re)**: BEAST 2 BDSKY
5. **Recombination Rate**: BEAST 2 CoalRe
6. **Genetic Diversity (π)**: Tajima
7. **Genetic Distance (GD)**
8. **Mutation Burden (MB)**
9. **GC Content**

*(Full methodological details are available directly in the Web Dashboard UI).*

## Acknowledgments & Citations

This tool acts as a powerful orchestrator for several state-of-the-art phylogenetic software packages. If you use Epi-GVI in your research, please ensure you cite both this repository and the incredible teams behind the underlying engines:

* **IQ-TREE 2** (Phylogenetic Inference & ModelFinder): 
  * *Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R (2020). IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era. Mol. Biol. Evol., 37:1530-1534.*
* **LSD2** (Least-Squares Dating / Evolutionary Rates): 
  * *To TH, Jung M, Ly-Trong N, et al. (2016). Fast dating using least-squares. Systematic Biology, 65(1):82-97.*
* **BEAST 2** (Bayesian MCMC Framework): 
  * *Bouckaert R, Vaughan TG, Barido-Sottani J, et al. (2019). BEAST 2.5: An advanced software platform for Bayesian evolutionary analysis. PLoS computational biology, 15(4), e1006650.*
* **BDSKY** (Birth-Death Skyline / Effective Reproductive Number):
  * *Stadler T, Kühnert D, Bonhoeffer S, Drummond AJ (2013). Estimating the past dynamics of macroevolution and epidemiology from molecular sequences. PNAS, 110(1):228-233.*
* **CoalRe** (Recombination Rates via Ancestral Recombination Graphs):
  * *Vaughan TG, Welch D, Drummond AJ, Standfield PJ, Hassler GW, Schwartz R (2014). Inferring ancestral recombination graphs from bacterial genomic data. Genetics, 198(1):257-270.*
