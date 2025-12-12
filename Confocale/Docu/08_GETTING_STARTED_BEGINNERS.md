# Getting Started Guide - MSD Analysis for Beginners

## Welcome!

This guide will help you get started with the MSD (Mean Squared Displacement) analysis toolkit, even if you've never used Python or command-line tools before.

---

## What This Toolkit Does

This software analyzes particle tracking data to understand how particles move. It can:

1. **Read your trajectory data** from CSV files
2. **Calculate movement statistics** (MSD - how far particles spread over time)
3. **Create plots** showing particle movement patterns
4. **Extract physical parameters** like diffusion coefficients
5. **Generate reports** summarizing the analysis

## What You Need

### 1. Your Data
- CSV files containing particle trajectory data
- Each file should have columns for: Track ID, X Position, Y Position, Time
- Example: Data exported from ImageJ TrackMate or similar tracking software

### 2. Software Requirements
- **Python 3.7+** (free programming language)
- **Required Python packages**: pandas, numpy, scipy, matplotlib

### 3. Basic Computer Skills
- Opening Command Prompt/PowerShell (we'll show you how)
- Navigating to folders
- Running simple commands

---

## Step-by-Step Setup

### Step 1: Install Python

1. **Download Python** from https://python.org/downloads/
2. **During installation**:
   - ✅ Check "Add Python to PATH" 
   - ✅ Check "Install pip"
   - Choose "Install Now"

3. **Verify installation**:
   - Press `Windows + R`, type `cmd`, press Enter
   - Type: `python --version`
   - You should see something like: `Python 3.11.5`

### Step 2: Install Required Packages

In the same Command Prompt window, type these commands one by one:

```cmd
pip install pandas
pip install numpy
pip install scipy
pip install matplotlib
```

Wait for each command to complete before typing the next one.

### Step 3: Download the Analysis Scripts

If you haven't already, make sure you have all the Python files in a folder on your computer.

### Step 4: Prepare Your Data

1. **Create the data folders**:
   - In the same folder as the Python scripts, create a folder called `Data`
   - Inside `Data`, create two folders:
     - `31_10_no_anomalous` (for normal diffusion data)
     - `17_11_anomalous` (for anomalous diffusion data)

2. **Put your CSV files** in the appropriate folder based on what type of motion you expect

---

## Your First Analysis

### Option 1: Complete Automatic Analysis (Recommended for Beginners)

1. **Open Command Prompt in the correct folder**:
   - Navigate to the folder containing the Python scripts
   - Hold `Shift` and right-click in an empty area
   - Select "Open PowerShell window here" or "Open Command Prompt here"

2. **Run the automated analysis**:
   
   **For normal diffusion data**:
   ```powershell
   python batch_analysis_no_anomalous.py
   ```
   
   **For anomalous diffusion data**:
   ```powershell
   python batch_analysis_anomalous.py
   ```

3. **Wait for completion** - this may take several minutes depending on your data size

4. **Check your results**:
   - New folders will be created with plots and analysis results
   - A summary report will be created in the `Docu` folder
   - Open the summary report (`.md` file) with any text editor

### Option 2: Analyze a Single File (For Learning)

1. **Create a simple plot** of your data:
   ```powershell
   python eamsd_plot.py "Data/31_10_no_anomalous/your_file.csv"
   ```
   Replace `your_file.csv` with the actual name of your file.

2. **Look at the results**:
   - A plot file (`.svg`) will be created
   - Open it with any web browser or image viewer

3. **Try fitting** to extract diffusion coefficient:
   ```powershell
   python msd_fitting.py "Data/31_10_no_anomalous/your_file.csv"
   ```

---

## Understanding Your Results

### Plots
- **MSD plots** show how particle spreading increases with time
- **Linear plots** indicate normal diffusion
- **Curved plots** may indicate anomalous diffusion or drift

### Key Numbers to Look For

#### Diffusion Coefficient (D)
- **Units**: μm²/s (micrometers squared per second)
- **Typical values for particles in glycerol**: 0.001 - 0.01 μm²/s
- **Meaning**: How fast particles spread due to random motion

#### R² Value
- **Range**: 0 to 1 (higher is better)
- **Good fit**: R² > 0.95
- **Meaning**: How well the mathematical model fits your data

#### Anomalous Exponent (α) - for anomalous data
- **α < 1**: Subdiffusion (particles are confined or hindered)
- **α = 1**: Normal diffusion (standard random walk)
- **α > 1**: Superdiffusion (particles have active transport)

### Summary Reports
- Located in the `Docu` folder
- Contains tables with all fitted parameters
- Easy to copy into Excel or other programs

---

## Common Beginner Mistakes & Solutions

### "Command not found" Error
**Problem**: `'python' is not recognized as an internal or external command`

**Solution**: 
1. Reinstall Python and make sure to check "Add Python to PATH"
2. Try using `py` instead of `python`:
   ```cmd
   py batch_analysis_no_anomalous.py
   ```

### "No module named pandas" Error
**Problem**: Required packages not installed

**Solution**: Install packages:
```cmd
pip install pandas numpy scipy matplotlib
```

### "No such file or directory" Error
**Problem**: File path is wrong

**Solutions**:
1. Make sure you're in the correct folder (where the Python scripts are)
2. Use the full path to your file:
   ```powershell
   python eamsd_plot.py "C:\Users\YourName\Documents\Data\your_file.csv"
   ```
3. Check that your file actually exists and is spelled correctly

### "No tracks found" Error
**Problem**: Your CSV file format isn't recognized

**Solution**: 
1. Open your CSV file in Excel or Notepad
2. Make sure it has columns named something like:
   - Track ID (or similar)
   - X position (or X, Position X, etc.)
   - Y position (or Y, Position Y, etc.)
   - Time (or T, Position T, etc.)
3. Check the [Data Requirements section](00_COMPLETE_USER_GUIDE.md#data-requirements) for supported column names

---

## What Each Folder Contains

After running the analysis, you'll see these new folders:

### `eamsd_plots/`
- Ensemble-averaged MSD plots
- Shows average behavior across all particles
- Linear axes plots

### `tamsd_plots/`
- Time-averaged MSD plots
- Shows behavior of individual particle tracks
- Good for understanding single-particle dynamics

### `linear_fits/` or `anomalous_analysis/linear_fits/`
- Results from fitting mathematical models
- Contains both plots (`.svg`) and numerical results (`.txt`)
- Diffusion coefficients and fit quality metrics

### `nonlinear_fits/` or `anomalous_analysis/nonlinear_fits/`
- Results from advanced models including drift/transport
- Contains drift velocities and anomalous exponents

### `Docu/`
- Summary reports in Markdown format (`.md`)
- Open with any text editor, Word, or web browser
- Tables with all fitted parameters organized by file

---

## Next Steps

### Once You're Comfortable with Basics

1. **Try different analysis parameters**:
   ```powershell
   # Use only first 30% of data for fitting
   python msd_fitting.py "your_file.csv" --manual-fraction 0.3
   
   # Analyze specific particle track
   python tamsd_plot.py "your_file.csv" --track-id 5
   ```

2. **Compare multiple experiments**:
   ```powershell
   python compare_msd.py
   ```
   Follow the interactive prompts to compare different datasets.

3. **Read the detailed documentation**:
   - `00_COMPLETE_USER_GUIDE.md` - Comprehensive user guide
   - `07_SCRIPT_REFERENCE.md` - Detailed technical documentation

### For Advanced Users

1. **Customize the batch analysis** by editing the scripts
2. **Write your own analysis scripts** using the provided modules
3. **Integrate with other analysis pipelines**

---

## Getting Help

### Built-in Help
Most scripts have help information:
```powershell
python eamsd_plot.py --help
python msd_fitting.py --help
```

### Check Your Work
Use the validation script:
```powershell
python test_msd_simulation.py
```
This creates simulated data with known parameters and checks if the analysis recovers them correctly.

### Documentation Files
1. **This file** - Basic getting started
2. **Complete User Guide** - Comprehensive documentation with examples
3. **Script Reference** - Technical details for each script

### Common Issues
See the [Troubleshooting section](00_COMPLETE_USER_GUIDE.md#troubleshooting) in the Complete User Guide.

---

## Example Workflow for Complete Beginners

Here's a complete example of analyzing your first dataset:

### 1. Setup (One-time only)
```cmd
# Install Python packages (if not done already)
pip install pandas numpy scipy matplotlib
```

### 2. Prepare Data
- Put your CSV file in `Data/31_10_no_anomalous/experiment1.csv`

### 3. Run Analysis
```powershell
# Open PowerShell in the script folder
python batch_analysis_no_anomalous.py
```

### 4. Check Results
- Look for new folders: `eamsd_plots`, `tamsd_plots`, `linear_fits`, `nonlinear_fits`
- Open `Docu/05_BATCH_ANALYSIS_LOG.md` to see summary table
- Open some `.svg` plot files in your web browser

### 5. Interpret Results
- Find your file in the summary table
- Note the Diffusion Coefficient (D) values
- Check R² values (should be > 0.9 for good fits)
- Look at plots to see if they make sense

### 6. Try Individual Analysis
```powershell
# Create a simple plot
python eamsd_plot.py "Data/31_10_no_anomalous/experiment1.csv"

# Try fitting with different fraction
python msd_fitting.py "Data/31_10_no_anomalous/experiment1.csv" --manual-fraction 0.25
```

---

## Success! You're Now Ready to Analyze Your Data

This toolkit provides powerful analysis capabilities for particle tracking data. Start with the automated batch analysis, then explore individual scripts as you become more comfortable.

Remember: 
- **Start simple** with batch analysis
- **Check the documentation** when you need more details
- **Experiment with parameters** once you understand the basics
- **Look at the plots** - they often tell you more than the numbers alone

Good luck with your analysis!