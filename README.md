<img width="560" height="420" alt="untitled1" src="https://github.com/user-attachments/assets/f85aee53-ff61-4528-912c-665e174f7562" />
<img width="560" height="420" alt="untitled" src="https://github.com/user-attachments/assets/ad68a715-5a97-4f64-b779-8509825dab59" />

# Adaptive Augmented Unscented Kalman Smoother for ECG Denoising

MATLAB implementation of an **adaptive augmented Unscented Kalman Filter / Smoother (AUKF/AUKS)** for ECG denoising with **non-additive noise**.

This project is based on the modeling framework proposed in:

**Efficient Bayesian ECG denoising using adaptive covariance estimation and nonlinear Kalman Filtering**  
Hamed Danandeh Hesar, Amin Danandeh Hesar  
*Computers and Electrical Engineering*, 2024  
DOI: [10.1016/j.compeleceng.2024.109869](https://doi.org/10.1016/j.compeleceng.2024.109869)

---

## Overview

ECG signals are often contaminated by:

- baseline wander
- muscle artifacts
- sensor noise
- power-line interference

This project implements a **Bayesian nonlinear state-space denoising framework** based on:

- phase-domain ECG modeling
- Gaussian mixture approximation of ECG morphology
- adaptive covariance estimation
- Unscented Kalman filtering
- Unscented Kalman smoothing

The code estimates a clean ECG signal from noisy measurements while continuously adapting the process and measurement noise covariances.

---

## Main Features

- **Adaptive Unscented Kalman Filter (AUKF)**
- **Adaptive Unscented Kalman Smoother (AUKS)**
- **ECG phase modeling** using R-peak locations
- **Gaussian mixture model** for ECG morphology
- **Online covariance adaptation** for \(Q\) and \(R\)
- **Backward smoothing** for improved denoising
- **SNR improvement evaluation**
- MATLAB implementation suitable for biomedical signal processing research

---

## Method Summary

The algorithm works in the following steps:

### 1. Load ECG Data
The code loads a `.mat` file containing:

- `x` → ECG signal matrix
- `fs` → sampling frequency

The first channel of `x` is used as the ECG signal.

---

### 2. Add Synthetic Noise
White Gaussian noise is added to the ECG signal at a selected SNR level:

```matlab
SNR = 6;
x_noisy(1,:) = awgn(x(1,:), SNR, 'measured');
```

This is used to test the denoising performance.

---

### 3. Detect R-Peaks
The code uses the **Pan–Tompkins algorithm**:

```matlab
[qrs_positions] = pantompkins_qrs(x_noisy(1,:), fs);
```

These R-peaks are used to define cardiac cycles.

---

### 4. Compute Linear ECG Phase
A phase signal is constructed from the RR intervals:

- phase is reset at each R-peak
- phase evolves linearly between peaks
- phase values are wrapped to \([-\pi, \pi]\)

This phase representation is useful for periodic ECG modeling.

---

### 5. Extract Mean ECG Morphology
The noisy ECG is binned by phase and the mean morphology is estimated:

```matlab
[ECGsd, ECGmean, meanphase] = ecgsd_extractor_ver1(...)
```

Then the mean ECG waveform is smoothed using wavelet denoising:

```matlab
ECGmean = wdenoise(...)
```

---

### 6. Fit Gaussian Mixture ECG Model
The mean ECG shape is approximated using a sum of Gaussian kernels:

```math
ECG(\theta) = \sum_{i=1}^{N} a_i \exp\left(-\frac{(\theta-\theta_i)^2}{2b_i^2}\right)
```

where:

- \(a_i\) = amplitude
- \(b_i\) = width
- \(theta_i\) = center phase

The parameters are estimated using **Particle Swarm Optimization (PSO)**.

Only the strongest Gaussian components are kept.

---

### 7. Estimate Instantaneous Angular Frequency
The angular frequency is derived from the RR intervals and stored as an additional measurement channel.

The measurement vector is:

```matlab
y = [x_noisy(1,:); x_noisy(2,:); x_noisy(3,:)];
```

where:

- channel 1 → noisy ECG
- channel 2 → phase
- channel 3 → instantaneous angular frequency

---

### 8. Adaptive Unscented Kalman Filtering
The filter uses an augmented state vector containing:

- ECG amplitude
- phase
- angular frequency
- Gaussian parameters
- process-noise augmentation terms

Sigma points are generated using the Unscented Transform.

The code uses a special configuration:

- `alpha = 1`
- `kappa = 0`
- `beta = 2`

This corresponds to a UKF/UKS-style implementation tailored for the augmented nonlinear ECG model.

---

### 9. Adaptive Covariance Update
The covariances are updated online:

- measurement covariance \(R\)
- process covariance \(Q\)

A sliding memory of innovations is used with a forgetting factor:

```matlab
forgetting_factor = 0.99;
window_size = round(fs/3);
```

This allows the filter to adapt to signal variability and changing noise characteristics.

---

### 10. Backward Smoothing
After the forward filter pass, the algorithm runs a **backward smoother** to refine the estimate:

- `AUKF` → forward estimate
- `AUKS` → smoothed estimate

The smoother usually produces a cleaner and more stable reconstruction.

---

## Input

The script expects a MATLAB `.mat` file with:

- `x` → ECG signal matrix
- `fs` → sampling frequency

Example:

```matlab
x = [ECG_signal];
fs = 360;
```

You will be prompted to select the file using a file dialog.

---

## Output

The code produces:

### 1. Adaptive UKF denoised ECG
Stored in:

```matlab
Xukf_update
```

### 2. Adaptive UKS smoothed ECG
Stored in:

```matlab
Xuks
```

### 3. SNR metrics
Computed as:

```matlab
AUKF_SNR
AUKS_SNR
```

These values indicate denoising performance.

---

## Figures Generated

### Figure 1
Noisy ECG signal with detected R-peaks.

### Figure 2
Comparison of:

- Original ECG
- Noisy ECG
- Adaptive UKF output

### Figure 3
Comparison of:

- Original ECG
- Noisy ECG
- Adaptive UKS output

---

## Key Parameters

Important tunable parameters in the script:

```matlab
SNR               % Added noise level
ecg_bins          % Number of phase bins
MaxNumGaussian    % Maximum number of Gaussian kernels retained
window_size       % Innovation memory window
forgetting_factor  % Adaptive covariance memory factor
alpha             % UT scaling parameter
beta              % UT prior knowledge parameter
kappa             % UT spread parameter
```

---

## Required MATLAB Functions

The following helper function must be available in the MATLAB path:

- `pantompkins_qrs.m`

The code also uses MATLAB built-in or toolbox functions such as:

- `awgn`
- `wdenoise`
- `particleswarm`

---

## MATLAB Toolboxes Required

- Signal Processing Toolbox
- Wavelet Toolbox
- Global Optimization Toolbox
- Communications Toolbox or equivalent support for `awgn`

---

## Applications

This implementation is useful for:

- ECG denoising research
- biomedical signal processing
- model-based cardiac signal analysis
- wearable ECG monitoring
- preprocessing for arrhythmia detection
- academic comparison of Kalman-based denoising methods

---

## Reference

**Primary paper:**

Hamed Danandeh Hesar, Amin Danandeh Hesar,  
*Efficient Bayesian ECG denoising using adaptive covariance estimation and nonlinear Kalman Filtering*,  
*Computers and Electrical Engineering*, 2024.  
DOI: [10.1016/j.compeleceng.2024.109869](https://doi.org/10.1016/j.compeleceng.2024.109869)

**Related paper:**

Hamed Danandeh Hesar, Amin Danandeh Hesar,  
*Adaptive augmented cubature Kalman filter/smoother for ECG denoising*,  
*Biomedical Engineering Letters*, 2024.  
DOI: [10.1007/s13534-024-00362-7](https://doi.org/10.1007/s13534-024-00362-7)

---

