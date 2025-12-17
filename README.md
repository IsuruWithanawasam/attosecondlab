# Attosecond Pulse Generation and Characterization

This repository contains Python scripts and analysis tools for investigating High-Order Harmonic Generation (HHG) and the characterization of attosecond pulses generated at the ELI ALPS facility.

## Overview

The code supports experimental work on attosecond pulse generation through HHG in noble gases, with particular focus on:
- IR driving laser beam characterization
- XUV radiation spectral analysis
- Spatial beam profile measurements
- Energy calibration and spectral characterization using filter transmission

## Repository Structure

```
attosecondlab/
├── beam_prof_analysis.py      # IR beam profile analysis
├── spec_analysis.py            # XUV spectral analysis and calibration
├── Report_Applications_of_attosecond_pulses.pdf  # Full experimental report
└── README.md
```

## Scripts

### 1. `beam_prof_analysis.py`

Analyzes the infrared (IR) driving laser beam profile recorded at the Holey Splitting Mirror (HSM).

**Features:**
- Gaussian profile fitting for horizontal and vertical beam cross-sections
- Background noise subtraction
- Beam waist calculation at the focal point
- Peak intensity estimation
- Energy loss correction due to beam truncation at HSM
- Ponderomotive energy (Up) calculation
- Cutoff energy estimation

**Key Parameters:**
- Wavelength (λ₀): 1030 nm
- Pulse energy: 1 mJ
- Pulse duration: 40 fs (FWHM)
- Focal length: 90 cm
- HSM hole diameter: 6 mm (for calibration)

**Outputs:**
- Beam radius at input (win) and focal point (w)
- Peak intensity (I₀): ~8.47 × 10¹⁴ W/cm²
- Ponderomotive energy: ~83.84 eV
- Theoretical cutoff energy: ~284.04 eV (for Argon)

**Usage:**
```python
python beam_prof_analysis.py
```

**Required Input:**
- `beam.png` - Beam profile image from HSM

### 2. `spec_analysis.py`

Performs comprehensive analysis of XUV spectra obtained from the Flat-Field Spectrometer (FFS).

**Features:**
- Background subtraction from dark regions
- Spectral extraction from 2D FFS images
- Transmission curve calculation using Tellurium (Te) filter
- Energy axis calibration using filter absorption edge
- Harmonic peak identification
- Spatial beam profile analysis
- XUV beam divergence calculation

**Key Analysis Steps:**
1. **Image Processing**: Background subtraction and noise reduction
2. **Spectral Extraction**: Vertical integration to obtain 1D spectrum
3. **Transmission Analysis**: Comparison of filtered vs. unfiltered spectra
4. **Energy Calibration**: 
   - Uses Te filter absorption edge at 43.1 eV as reference
   - Fits calibration curve: E = a/x + b
   - Harmonic spacing: 2.408 eV (2 × photon energy)
5. **Spatial Analysis**: Horizontal integration for beam profile

**Outputs:**
- Calibrated XUV spectrum
- Macroscopic cutoff energy: ~67.44 eV
- XUV beam divergence: 0.0012 rad
- Divergence ratio (θ_XUV/θ_IR): 0.126

**Usage:**
```python
python spec_analysis.py
```

**Required Inputs:**
- FFS image without filter: `task13-15_XUV_flatfield_NO_filter_*.png`
- FFS image with Te filter: `task13-15_XUV_flatfield_Te100nm_filter_*.png`
- Filter transmission data: `filter.txt`

## Theoretical Background

### High-Order Harmonic Generation (HHG)

HHG is a nonlinear optical process where intense IR laser pulses interact with atomic gases to produce coherent XUV radiation. The process is described by the **three-step model**:

1. **Ionization**: Strong laser field tunnel-ionizes valence electrons
2. **Acceleration**: Free electrons gain kinetic energy in the oscillating field
3. **Recombination**: Electrons return and recombine, emitting XUV photons

### Key Equations

**Cutoff Energy:**
```
E_cutoff = Ip + 3.2 × Up
```
where:
- Ip = ionization potential (15.76 eV for Argon)
- Up = ponderomotive potential

**Ponderomotive Potential:**
```
Up [eV] = 9.33 × 10⁻¹⁴ × IL [W/cm²] × (λ₀ [μm])²
```

**Focal Spot Size:**
```
w = (λ₀/π) × (f/win)
```

**Peak Intensity:**
```
I₀ = (2εn)/(τ × π × w²) × √(4ln2/π)
```

## Key Findings

### Microscopic vs. Macroscopic Cutoff
- **Microscopic cutoff** (single-atom, theoretical): 284.04 eV
- **Macroscopic cutoff** (experimental): 67.44 eV
- **Discrepancy**: Primarily due to **phase-matching effects** in the macroscopic gas medium

### Beam Divergence
The XUV beam shows significantly tighter collimation than the IR driving beam:
- θ_XUV = 0.0012 rad
- θ_IR = 0.0091 rad
- Ratio = 0.126

This confirms that XUV radiation inherits favorable spatial coherence properties from the driving laser.

## Dependencies

```bash
pip install numpy matplotlib scipy pillow
```

**Required packages:**
- `numpy` - Numerical computations
- `matplotlib` - Data visualization
- `scipy` - Signal processing and curve fitting
- `PIL` (Pillow) - Image loading and processing

## Installation

```bash
# Clone the repository
git clone https://github.com/IsuruWithanawasam/attosecondlab.git
cd attosecondlab

# Install dependencies
pip install -r requirements.txt
```

## Experimental Setup

**Location:** ELI ALPS HR Condensed Beamline

**Laser Parameters:**
- Source: High-average-power femtosecond laser
- Wavelength: 1030 nm
- Pulse energy: 1 mJ
- Pulse duration: 40 fs (FWHM)
- Repetition rate: High-average-power

**Target Gas:** Argon (Ar)
- Ionization potential: 15.76 eV

**Diagnostics:**
- Flat-Field Spectrometer (FFS) for XUV spectral measurements
- MCP detector with Te filter option
- Beam profiler at Holey Splitting Mirror

## Data Analysis Workflow

```
1. Beam Profile Analysis
   └── beam_prof_analysis.py
       ├── Load HSM beam image
       ├── Gaussian fitting (H & V)
       ├── Calculate focal parameters
       └── Estimate intensity & cutoff

2. Spectral Analysis
   └── spec_analysis.py
       ├── Load FFS images (with/without filter)
       ├── Background subtraction
       ├── Extract spectra
       ├── Calculate transmission
       ├── Energy calibration
       └── Spatial profile analysis
```

## Results Visualization

The scripts generate multiple visualization plots:

1. **Beam Profile**: 2D image with intensity cross-sections
2. **Gaussian Fits**: Horizontal and vertical intensity profiles
3. **XUV Spectra**: Raw and calibrated photon energy spectra
4. **Transmission Curve**: Te filter transmission vs. pixel/energy
5. **Calibration Curve**: Pixel-to-energy mapping with error bars
6. **Spatial Profile**: Vertical beam distribution

## Configuration

Key parameters can be modified in the scripts:

**beam_prof_analysis.py:**
```python
IMAGE_FILE = 'beam.png'
HLINE = 633          # Horizontal profile line
VLINE = 685          # Vertical profile line
d = 6e-3             # HSM hole diameter
lmd = 1030e-9        # Wavelength
f = 90e-2            # Focal length
```

**spec_analysis.py:**
```python
FILE_NO_FILTER = 'task13-15_XUV_flatfield_NO_filter_*.png'
FILE_FILTER = 'task13-15_XUV_flatfield_Te100nm_filter_*.png'
DARK_REGION_COORDS = [800, 1200, 0, 500]
SPECTRUM_PIXEL_LIMITS = [500, 1500]
```

## Calibration Notes

### Energy Calibration Methodology
The spectral energy axis is calibrated using the characteristic absorption edge of the Tellurium filter:
- **Reference energy**: 43.1 eV (Te M-edge)
- **Calibration model**: E(x) = a/x + b
- **Harmonic spacing**: 2.408 eV (twice the fundamental photon energy at 1030 nm)
- **Fitting quality**: R² = 1.00

Error sources considered:
- Pixel localization uncertainty: ±0.3 pixels
- Filter edge position uncertainty
- Peak identification accuracy

## Troubleshooting

**Common Issues:**

1. **File not found errors**
   - Ensure input images are in the same directory as the scripts
   - Check filenames match exactly

2. **Peak detection fails**
   - Adjust `height` parameter in `find_peaks()`
   - Check background subtraction coordinates

3. **Poor Gaussian fits**
   - Verify exclusion ranges (H_LIM, V_LIM)
   - Check for saturation or clipping in images

4. **Calibration errors**
   - Verify filter transmission data file format
   - Check harmonic spacing matches laser parameters

## References

1. Chang, Z. (2011). *Fundamentals of Attosecond Optics*. CRC Press.
2. Major, B. & Tímár-Grósz, T. (2025). *Applications of Attosecond Pulses: Practice Syllabus*. ELI ALPS.

## Citation

If you use this code in your research, please cite:

```bibtex
@misc{attosecondlab2025,
  author = {Mawella Withanawasam, Isuru Chanilka},
  title = {Attosecond Pulse Generation and Characterization Analysis Tools},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/IsuruWithanawasam/attosecondlab}
}
```

## Author

**Isuru Chanilka Mawella Withanawasam**
- LAScALA Erasmus Mundus Master Program
- ELI ALPS Facility
- Experiment Date: October 9, 2025

## License

This project is available for academic and research purposes. Please contact the author for commercial use.

## Acknowledgments

- ELI ALPS HR Condensed Beamline team
- LAScALA Erasmus Mundus Master Program
- Supervisors and laboratory staff at ELI ALPS

---

**Last Updated:** October 30, 2025

For questions or issues, please open an issue on GitHub or contact the author directly.
