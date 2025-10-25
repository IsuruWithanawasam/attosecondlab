# -*- coding: utf-8 -*-
"""
Gaussian Beam Profile Analysis
Analyzes horizontal and vertical intensity profiles of a beam image
Created on Fri Oct 10 22:20:00 2025
@author: Asus
"""

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, median_filter
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

# ============================================================================
# CONFIGURATION
# ============================================================================

IMAGE_FILE = 'beam.png'
HLINE = 633  # Horizontal line for profile extraction
VLINE = 685  # Vertical line for profile extraction

# Horizontal profile exclusion range
H_LIM_LOW = 615
H_LIM_HIGH = 770

# Vertical profile exclusion range
V_LIM_LOW = 565
V_LIM_HIGH = 710

d=6e-3
lmd=1030e-9
f=90e-2
e_n=1e-3
tau=40e-15

# ============================================================================
# FUNCTIONS
# ============================================================================

def enReduction(I0_initial,w_in):
    # --- Data extracted from Figure 4(b) ---
    beam_sizes_mm = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]) # x-axis
    energy_in_hole_percent = np.array([100, 98, 85, 58, 30, 15, 8, 5, 4, 3, 2]) # y-axis

    # --- Create an interpolation function from the data ---
    # This function can now estimate the % loss for ANY beam size in the range
    get_loss_percentage = interp1d(beam_sizes_mm, energy_in_hole_percent, kind='cubic')

    # --- Your previously calculated value ---
    w_in_calculated_mm = w_in*1e3 # Example value

    # --- Step 1: Use the function to find the loss ---
    percent_energy_lost_accurate = get_loss_percentage(w_in_calculated_mm)

    # --- The rest of the calculation is the same ---
    original_pulse_energy_J = 1e-3
    energy_fraction_remaining = (100.0 - percent_energy_lost_accurate) / 100.0
    reduced_pulse_energy_J = original_pulse_energy_J * energy_fraction_remaining

    # Recalculate intensity as before...
    I0_reduced = I0_initial * energy_fraction_remaining

    print(f"For a beam size of {w_in_calculated_mm} mm:\nAccurate percentage of energy lost: {percent_energy_lost_accurate:.2f}%\nCorrected peak intensity would be scaled by a factor of {energy_fraction_remaining:.2f}")
    
    return I0_reduced

def gaussian_with_bg(x, A, x0, w, B):
    """
    Gaussian function with background offset.
    
    Parameters:
        x: position
        A: amplitude
        x0: center position
        w: width parameter
        B: background offset
    """
    return A * np.exp(-2 * (x - x0)**2 / w**2) + B


def load_and_preprocess_image(filename):
    """Load image and apply filtering."""
    img = Image.open(filename)
    print(f"Image mode: {img.mode}")
    
    # Apply median filter followed by Gaussian blur
    img_filtered = gaussian_filter(median_filter(img, size=3), sigma=1.2)
    img_matrix = np.array(img_filtered)
    
    return img_matrix


def extract_and_plot_profiles(img_matrix, hline, vline):
    """Extract horizontal and vertical profiles and display."""
    height, width = img_matrix.shape
    print(f"Image dimensions: {height} x {width}")
    
    # Display image with profile lines
    plt.figure()
    plt.imshow(img_matrix)
    plt.xlabel('Pixel (X)')
    plt.ylabel('Pixel (Y)')
    plt.axhline(hline, color='red', label='Horizontal profile line')
    plt.axvline(vline, color='green', label='Vertical profile line')
    plt.legend()
    plt.show()
    
    # Extract profiles
    h_profile = img_matrix[hline, :]
    v_profile = img_matrix[:, vline]
    
    return h_profile, v_profile, width, height


def fit_and_analyze_profile(profile, pixel_indices, lim_low, lim_high, direction):
    """
    Fit Gaussian to profile data and display results.
    
    Parameters:
        profile: intensity profile data
        pixel_indices: pixel position indices
        lim_low, lim_high: range to exclude from fit
        direction: 'Horizontal' or 'Vertical' (for labeling)
    
    Returns:
        tuple: fitted parameters (A, x0, w, B)
    """
    # Remove excluded region
    pixel_indices_trimmed = np.delete(pixel_indices, np.arange(lim_low, lim_high))
    profile_trimmed = np.delete(profile, np.arange(lim_low, lim_high))
    
    # Initial parameter guesses
    p0 = [
        np.max(profile_trimmed) - np.min(profile_trimmed),  # Amplitude
        np.mean(pixel_indices_trimmed),                      # Center
        100,                                                 # Width
        np.min(profile_trimmed)                              # Background
    ]
    
    # Fit curve
    popt, pcov = curve_fit(gaussian_with_bg, pixel_indices_trimmed, profile_trimmed, p0=p0)
    A_fit, x0_fit, w_fit, B_fit = popt
    
    # Plot results
    plt.figure()
    plt.plot(pixel_indices_trimmed, profile_trimmed, '.', label='Data')
    plt.plot(pixel_indices, gaussian_with_bg(pixel_indices, *popt), 'r-', label='Gaussian fit')
    plt.xlabel("Pixel position")
    plt.ylabel("Intensity")
    plt.title(f"{direction} Profile - Gaussian Fit")
    plt.legend()
    plt.show()
    
    # Print results
    print(f"\n{direction} Profile Results:")
    print(f"  Amplitude (A): {A_fit:.4f}")
    print(f"  Center (x0): {x0_fit:.4f}")
    print(f"  Width (w): {w_fit:.4f}")
    print(f"  Background (B): {B_fit:.4f}")
    
    return A_fit, x0_fit, w_fit, B_fit


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

# Load and preprocess image
img_matrix = load_and_preprocess_image(IMAGE_FILE)
height, width = img_matrix.shape

# Extract and display profiles
h_profile, v_profile, width, height = extract_and_plot_profiles(img_matrix, HLINE, VLINE)

# Analyze horizontal profile
print("\n" + "="*60)
print("HORIZONTAL PROFILE ANALYSIS")
print("="*60)
pixH = np.arange(width)
A_h, x0_h, w_h, B_h = fit_and_analyze_profile(
    h_profile, pixH, H_LIM_LOW, H_LIM_HIGH, "Horizontal"
)

# Analyze vertical profile
print("\n" + "="*60)
print("VERTICAL PROFILE ANALYSIS")
print("="*60)
pixV = np.arange(height)
A_v, x0_v, w_v, B_v = fit_and_analyze_profile(
    v_profile, pixV, V_LIM_LOW, V_LIM_HIGH, "Vertical"
)


# ============================================================================
# COMBINED VISUALIZATION
# ============================================================================

print("\n" + "="*60)
print("GENERATING COMBINED VISUALIZATION")
print("="*60)

# Prepare data for subplots
pixH_trimmed = np.delete(np.arange(width), np.arange(H_LIM_LOW, H_LIM_HIGH))
h_profile_trimmed = np.delete(h_profile, np.arange(H_LIM_LOW, H_LIM_HIGH))
pixV_trimmed = np.delete(np.arange(height), np.arange(V_LIM_LOW, V_LIM_HIGH))
v_profile_trimmed = np.delete(v_profile, np.arange(V_LIM_LOW, V_LIM_HIGH))

# Create figure with custom layout
fig = plt.figure(figsize=(16, 12))

# Define grid
gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.35)

# Horizontal profile subplot (top-left, spans 2 columns)
ax_horiz = fig.add_subplot(gs[0, :2])
ax_horiz.plot(pixH_trimmed, h_profile_trimmed, 'g.', markersize=8, label='Data', alpha=0.7)
ax_horiz.plot(np.arange(width), gaussian_with_bg(np.arange(width), A_h, x0_h, w_h, B_h), 
              'r-', linewidth=3, label='Gaussian fit')
ax_horiz.set_ylabel("Intensity", fontsize=12, fontweight='bold')
ax_horiz.set_title("Horizontal Profile", fontsize=13, fontweight='bold')
ax_horiz.legend(loc='best', fontsize=11)
ax_horiz.grid(True, alpha=0.3, linestyle='--')
ax_horiz.set_xlim([width, 0])  # Match image x-axis range and direction
ax_horiz.tick_params(labelbottom=True)  # Hide x-axis labels for alignment
ax_horiz.invert_xaxis()

# Main image subplot (middle-left and middle-center, spans 2x2)
ax_main = fig.add_subplot(gs[1:3, :2])
ax_main.imshow(img_matrix)
ax_main.axhline(HLINE, color='red', linewidth=2.5, alpha=0.8, label='H-profile line')
ax_main.axvline(VLINE, color='blue', linewidth=2.5, alpha=0.8, label='V-profile line')
ax_main.set_xlabel("Pixel (X)", fontsize=12, fontweight='bold')
ax_main.set_ylabel("Pixel (Y)", fontsize=12, fontweight='bold')
ax_main.set_title("Beam Image with Profile Lines", fontsize=13, fontweight='bold')
ax_main.legend(loc='upper right', fontsize=11)

# Vertical profile subplot (right side, spans 2 rows)
ax_vert = fig.add_subplot(gs[1:3, 2])
ax_vert.plot(v_profile_trimmed, pixV_trimmed, 'b.', markersize=8, label='Data', alpha=0.7)
ax_vert.plot(gaussian_with_bg(pixV, A_v, x0_v, w_v, B_v), pixV, 'r-', 
             linewidth=3, label='Gaussian fit')
ax_vert.set_xlabel("Intensity", fontsize=12, fontweight='bold')
ax_vert.set_ylabel("Pixel (Y)", fontsize=12, fontweight='bold')
ax_vert.set_title("Vertical Profile", fontsize=13, fontweight='bold')
ax_vert.legend(loc='best', fontsize=11)
ax_vert.grid(True, alpha=0.3, linestyle='--')
ax_vert.invert_yaxis()
ax_vert.set_ylim([height, 0])  # Match image y-axis range and direction

# Hide top-right subplot
ax_empty = fig.add_subplot(gs[0, 2])
ax_empty.axis('off')

plt.show()

print("\n" + "="*60)
print("Beam Analysis")
print("="*60)

########################## FWHM calibration horizontal
w_h_m=d*w_h/(H_LIM_HIGH-H_LIM_LOW) # horizontal input beam width radius in meter
w_v_m=d*w_v/(V_LIM_HIGH-V_LIM_LOW)

w_f_h_m=lmd*f/(np.pi*w_h_m)                # horizontal focused beam width radius in meter
w_f_v_m=lmd*f/(np.pi*w_v_m)                # vertical focused beam width radius in meter

print(f"  1/e^2 radius horizontal (input): {w_h_m*1e3:.4f} mm")
print(f"  1/e^2 radius vertical (input): {w_v_m*1e3:.4f} mm")

print(f"  1/e^2 radius horizontal (focused): {w_f_h_m*1e3:.4f} mm")
print(f"  1/e^2 radius vertical (focused): {w_f_v_m*1e3:.4f} mm")


#I_0= (2*e_n/tau*(w_f_h_m*w_f_v_m)*np.pi)*np.sqrt(4*np.log(2)/np.pi)
I_0 = (2 * e_n / (tau * np.pi * w_f_h_m * w_f_v_m)) * np.sqrt(4 * np.log(2) / np.pi)
print(f"  Peak intensity (focused): {I_0*1e-4:.4f} W/cm2")

I_0_reduced=enReduction(I_0,(w_h_m+w_v_m)/2)
print(f"  Peak intensity (reduced): {I_0_reduced*1e-4:.4f} W/cm2")

print("\n" + "="*60)
print("Energy calculations")
print("="*60)

up_ev_nonr = 9.33e-14 * I_0 * 1e-4 * ((lmd * 1e6)**2)
up_ev_reduced = 9.33e-14 * I_0_reduced * 1e-4 * ((lmd * 1e6)**2)

print(f"  Pondemotive energy (Up): {up_ev_nonr:.4f} eV")
print(f"  Reduced Pondemotive energy (Up): {up_ev_reduced:.4f} eV")

I_p=15.76;# ev # for argon
E_cutoff_nr=I_p+up_ev_nonr*3.2
E_cutoff_reduced=I_p+up_ev_reduced*3.2


print(f"  Cutoff energy (Up): {E_cutoff_nr:.4f} eV")
print(f"  Reduced Cutoff energy (Up): {E_cutoff_reduced:.4f} eV")