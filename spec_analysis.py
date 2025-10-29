# -*- coding: utf-8 -*-
"""

Created on Fri Oct 10 15:19:35 2025
@author: Isuru Withanawasam
"""

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import math


# --- Input Files ---
FILE_NO_FILTER = 'task13-15_XUV_flatfield_NO_filter_Ar_80W_cellC_120mbar_abm11300_hm1140_selection7200_expt5ms_BB_s_pol_MCP6000000001.png'
FILE_FILTER = 'task13-15_XUV_flatfield_Te100nm_filter_Ar_80W_cellC_120mbar_abm11300_hm1140_selection7200_expt5ms_BB_s_pol_MCP600_Te100nmfilter0000001.png'
filterF='filter.txt'

DARK_REGION_COORDS = [800, 1200, 0, 500]
SPECTRUM_PIXEL_LIMITS = [500, 1500]


# Load image and convert to a numpy array for calculations
img_no_filter = Image.open(FILE_NO_FILTER)
img_matrix_no_filter = np.array(img_no_filter, dtype=float)

# Define the dark region from the coordinates
y_start, y_end, x_start, x_end = DARK_REGION_COORDS
dark_region_no_filter = img_matrix_no_filter[y_start:y_end, x_start:x_end]

# Calculate and subtract the background noise
background_avg_no_filter = np.mean(dark_region_no_filter)
print(f"  > Average background: {background_avg_no_filter:.2f}")
img_bg_subtracted_no_filter = img_matrix_no_filter - background_avg_no_filter

# Clip negative values to zero (intensity cannot be negative)
img_cleaned_no_filter = np.clip(img_bg_subtracted_no_filter, 0, None)

# Calculate the horizontal (spectral) profile by summing vertically
full_h_profile_no_filter = np.sum(img_cleaned_no_filter, axis=0)

# Apply the spectral pixel limits to get the final spectrum
#spec_start, spec_end = SPECTRUM_PIXEL_LIMITS
#S_wo = full_h_profile_no_filter[spec_start:spec_end]
S_wo = full_h_profile_no_filter

plt.figure()
plt.imshow(img_cleaned_no_filter)
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.show()



# Load image and convert to a numpy array
img_filter = Image.open(FILE_FILTER)
img_matrix_filter = np.array(img_filter, dtype=float)

# Define the dark region from the same coordinates
dark_region_filter = img_matrix_filter[y_start:y_end, x_start:x_end]

# Calculate and subtract the background noise
background_avg_filter = np.mean(dark_region_filter)
print(f"  > Average background: {background_avg_filter:.2f}")
img_bg_subtracted_filter = img_matrix_filter - background_avg_filter

# Clip negative values to zero
img_cleaned_filter = np.clip(img_bg_subtracted_filter, 0, None)

# Calculate the horizontal (spectral) profile
full_h_profile_filter = np.sum(img_cleaned_filter, axis=0)

# Apply the spectral pixel limits to get the final spectrum
#S_w = full_h_profile_filter[spec_start:spec_end]
S_w = full_h_profile_filter

# =============================================================================
# 5. CALCULATE TRANSMISSION AND PLOT RESULTS
# =============================================================================


# Calculate transmission
T = S_w / S_wo

# Create the plot with a primary and secondary y-axis
fig, ax1 = plt.subplots(figsize=(12, 7))

# Plot spectra on the primary (left) y-axis
color1 = 'tab:blue'
ax1.set_xlabel('Pixel Coordinate (within specified limits)')
ax1.set_ylabel('Integrated Intensity (counts)', color=color1)
ax1.plot(S_wo, color=color1, alpha=0.8, label='Spectrum without Filter ($S_{wo}$)')
ax1.plot(S_w, color='tab:orange', label='Spectrum with Filter ($S_w$)')
ax1.tick_params(axis='y', labelcolor=color1)
ax1.grid(True, linestyle='-', alpha=0.6)
ax1.set_ylim(bottom=0) # Ensure y-axis starts at 0

# Create the secondary (right) y-axis for transmission
ax2 = ax1.twinx()
color2 = 'tab:green'
ax2.set_ylabel('Transmission (T)', color=color2)
ax2.plot(T, color=color2, linestyle='-', label='Transmission (T = $S_w / S_{wo}$)')
ax2.tick_params(axis='y', labelcolor=color2)
ax2.set_ylim(0, max(1.0, np.percentile(T, 99.5))) # Set reasonable limits for transmission

# Add title and a combined legend
fig.suptitle('XUV Spectra and Experimental Transmission', fontsize=16)
fig.legend(loc='upper right', bbox_to_anchor=(0.9, 0.9))
fig.tight_layout(rect=[0, 0, 1, 0.96]) # Adjust layout to make room for title

plt.show()


min_T_pix=(np.where(T==min(T)))[0][0]
print('Position of the dip Transmission: ',min_T_pix)

##########################################
filter_d = np.loadtxt(filterF)
# Split into two arrays
filE = filter_d [:, 0]   # first column
T = filter_d [:, 1]   # second column

plt.plot(filE,T)
plt.xlabel('Photon energy(ev)')
plt.ylabel('Transmission')
plt.title('Tellurium Filter (density=6.24, Thickness=0.2 microns)')
plt.show()

# Select only the region
mask = (filE >= 30) & (filE <= 45)
E_region = filE[mask]
T_region = T[mask]

# Find index of minimum within this region
idx_min = np.argmin(T_region)

# Get corresponding x and y values
E_dip = E_region[idx_min]
T_dip = T_region[idx_min]


print('Energy of the dip filter: ',E_dip)
###########################

peaks,_ = find_peaks(S_wo,height=25000)
#print(peaks)
peaks=np.delete(peaks, [9,13])
plt.plot(range(len(S_wo)),S_wo)
plt.plot(peaks, S_wo[peaks], "x")
plt.xlabel('pixel')
plt.ylabel('Intensity count')
plt.show()

print('after filtering:',peaks)

################
c_peak= peaks[np.abs(peaks - min_T_pix).argmin()]

print('Closer peak position for the filter: ',c_peak)

# Define the energy separation between adjacent ODD harmonics
PHOTON_ENERGY_EV = 1.204
HARMONIC_SPACING_EV = 2 * PHOTON_ENERGY_EV
print(f"Energy spacing between harmonics: {HARMONIC_SPACING_EV:.3f} eV\n")

# =============================================================================
# 2. GENERATE CALIBRATION DATA POINTS
# =============================================================================
# Find which experimental peak is closest to the transmission dip's pixel.
# This will be our reliable anchor point for the energy scale.
peak_distances_from_dip = np.abs(peaks - c_peak)
anchor_peak_index = np.argmin(peak_distances_from_dip)
anchor_peak_pixel = peaks[anchor_peak_index]

print(f"The harmonic peak at pixel {anchor_peak_pixel} is the closest to the dip at pixel {min_T_pix}.")
print(f"This peak is our anchor and is assigned the energy {E_dip:.2f} eV.\n")

# Assign energies to all other peaks relative to the anchor peak
# The difference in index tells us how many harmonic steps away we are
harmonic_steps = np.arange(len(peaks)) - anchor_peak_index

# We assign energies based on these steps. Note the minus sign:
# smaller pixel numbers correspond to higher energies.
assigned_energies = E_dip - (harmonic_steps * HARMONIC_SPACING_EV)


##########################error calculation

pixel_spacing = np.mean(np.diff(peaks))        # average pixel distance between harmonic peaks
dE_dpix = HARMONIC_SPACING_EV / pixel_spacing  # energy change per pixel [eV/pixel]

# --- Step 2. Pixel uncertainties
delta_x_anchor = abs(c_peak - min_T_pix)       # pixel offset between chosen anchor and true dip
peak_pixel_localization_error = 0.3            # example uncertainty from peak fitting [pixels]

# --- Step 3. Convert to energy uncertainties
dE_anchor = delta_x_anchor * dE_dpix           # systematic calibration offset [eV]
dE_local = peak_pixel_localization_error * dE_dpix  # random localization uncertainty [eV]

# --- Step 4. Combine in quadrature
energy_errors = np.sqrt(dE_anchor**2 + dE_local**2)

# --- Step 5. Assign same uncertainty to all harmonics
assigned_energy_errors = np.full_like(assigned_energies, energy_errors)


"""
print("Generated (Pixel, Energy) pairs for calibration:")
for pix, nrg in zip(peaks, assigned_energies):
    print(f"  Pixel: {pix:4d} -> Energy: {nrg:5.2f} eV")
"""
# =============================================================================
# 3. FIT THE CALIBRATION CURVE
# =============================================================================
# Define the model function to be fitted: E = a/x + b
def calibration_model(x, a, b):
    return (a / x) + b

# Use curve_fit to find the best parameters 'a' and 'b'
popt, pcov = curve_fit(
    calibration_model,
    peaks,
    assigned_energies,
    sigma=assigned_energy_errors,
    absolute_sigma=True
)
a_fit, b_fit = popt
perr = np.sqrt(np.diag(pcov))  # standard deviation of each fitted parameter
a_err, b_err = perr

print(f"a = {a_fit:.5f} ± {a_err:.5f}")
print(f"b = {b_fit:.5f} ± {b_err:.5f}")

ss_res = np.sum((assigned_energies - calibration_model(peaks, a_fit, b_fit))**2)
ss_tot = np.sum((assigned_energies - np.mean(assigned_energies))**2)
R2 = 1 - ss_res/ss_tot
print(f"R² = {R2:.4f}")

# =============================================================================
# 4. VISUALIZE THE CALIBRATION FIT
# =============================================================================
# Create a smooth curve using the fitted parameters for visualization
pixel_range_for_plot = np.linspace(500, 1700, 1000)
fitted_curve = calibration_model(pixel_range_for_plot, a_fit, b_fit)

plt.figure(figsize=(10, 6))
plt.errorbar(
    peaks, 
    assigned_energies, 
    yerr=assigned_energy_errors,   # add the error bars
    fmt='ro',                      # red circles for the points
    ecolor='black',                # color of error bars
    capsize=3,                     # small horizontal caps on error bars
    label='Data Points (Harmonic Peaks)'
)
plt.plot(pixel_range_for_plot, fitted_curve, 'b-', label='Fitted Curve ($E = a/x + b$)')
plt.xlabel('Pixel Position')
plt.ylabel('Photon Energy (eV)')
plt.title('Spectrometer Calibration Curve Fit')

plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()


# =============================================================================
# 5. APPLY CALIBRATION TO THE FULL SPECTRUM
# =============================================================================
# Create a new energy axis for your entire spectrum using the fitted model
pixel_axis_full = np.arange(0,len(S_wo),1)
energy_axis_full = calibration_model(pixel_axis_full, a_fit, b_fit)



# Plot the final, calibrated spectrum

plt.plot(energy_axis_full[600:1700],S_wo[600:1700])
plt.xlabel('Photon Energy (eV)')
plt.ylabel('Intensity (counts)')
plt.title('Calibrated XUV Spectrum')
plt.grid(True)
# Invert x-axis so that higher energy is on the left (conventional for spectra)
#plt.gca().invert_xaxis()
plt.show()

#print(energy_axis_full)
#print(np.size(S_wo[1:len(S_wo)]))


########################################

# Calculate the horizontal (spectral) profile by summing vertically
full_v_profile_no_filter = np.sum(img_cleaned_no_filter, axis=1)
b_prof=full_v_profile_no_filter

v_lenght=np.size(b_prof)
vpix=np.arange(0,v_lenght)

phy_lenght=40 #mm

y=vpix*(phy_lenght/(v_lenght-1))

plt.plot(y,b_prof)
plt.xlabel('y (mm)')
plt.ylabel('Intensity (counts)')
plt.grid(True)
plt.show()



# Example: replace these with your actual arrays
x = np.array(y)  # x-values
I = np.array(b_prof)  # I-values

# Define the Gaussian function with offset
def gaussian_with_offset(x, A, x0, win, B):
    return A * np.exp(-2 * ((x - x0)**2) / win**2) + B

# Initial guess for parameters: [A, x0, win, B]
initial_guess = [max(I)-min(I), x[np.argmax(I)], (max(x)-min(x))/4, min(I)]

# Fit the curve
popt, pcov = curve_fit(gaussian_with_offset, x, I, p0=initial_guess)

# Extract fitted parameters
A_fit, x0_fit, win_fit, B_fit = popt
print("Fitted parameters:")
print(f"A = {A_fit}")
print(f"x0 = {x0_fit}")
print(f"win = {win_fit}")
print(f"B = {B_fit}")

# Compute fitted curve
I_fit = gaussian_with_offset(x, *popt)

# Plot original data and fitted curve
plt.plot(x, I, 'b', label='Data')
plt.plot(x, I_fit, 'r-', label='Fitted curve')
plt.xlabel('y (mm)')
plt.ylabel('Intensity (counts)')
#plt.title('Gaussian Fit with Offset')
plt.legend()
plt.show()

D=2.5



teta = win_fit * (1e-3) / D
theta_XUV = teta

print(f"Divergence of beam (rad) = {teta}")
print(f"Divergence of beam (XUV) (rad) = {theta_XUV}")

# Convert to degrees
teta_deg = teta * 180 / math.pi
theta_XUV_deg = theta_XUV * 180 / math.pi

#print(f"Divergence of beam (deg) = {teta_deg}")
#print(f"Divergence of beam (XUV) (deg) = {theta_XUV_deg}")

### from previous code
f = 90e-2  # m
w_h = 8.2799e-3  # m
w_v = 8.2082e-3  # m

theta_IR = w_v / f
print(f"Divergence of beam (IR) (rad) = {theta_IR}")

# Convert IR divergence to degrees
theta_IR_deg = theta_IR * 180 / math.pi
#print(f"Divergence of beam (IR) (deg) = {theta_IR_deg}")

# Ratio
ratio_XUV_IR = theta_XUV / theta_IR
print(f"Ratio (XUV/IR) = {ratio_XUV_IR}")


