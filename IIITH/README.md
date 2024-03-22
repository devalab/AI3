# AI3 - IIITH

## Prerequisites

Before setting up the Conda environment, ensure you have:

- Anaconda or Miniconda
- UCSF Chimera [Download](https://www.cgl.ucsf.edu/chimera/download.html)

## Setting Up the Conda Environment

### 1. Create Conda Environment

Use the `environment.yml` file to create a Conda environment with necessary dependencies:

```bash
cd IIITH
conda env create -f environment.yml
```

### 2. Activate the Environment

Activate the created environment:

```bash
conda activate aicube
```

## Running the Scripts

After setting up and activating the Conda environment:

### 1. Navigate to the Script Directory

```bash
cd scripts
```

### 2. Run the Script

```bash
python run_scripts.py
```
