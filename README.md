# Target Detection & Velocity Estimation using FMCW Radar

**Course:** EECG231 - Analysis of Continuous and Discrete-Time Signals (Fall 2025)
**Technology:** MATLAB
**Domain:** Automotive Radar, Signal Processing (DSP)

## Project Overview
[cite_start]This project simulates a **Frequency Modulated Continuous Wave (FMCW)** radar system operating in the **76.5 GHz automotive band**[cite: 324]. [cite_start]The system is designed to detect multiple moving targets and estimate their range and radial velocity with high precision using a complete signal processing chain[cite: 325].

[cite_start]The design implements a Linear Frequency Modulated (LFM) waveform and utilizes a **2D-FFT** architecture to generate Range-Doppler Maps[cite: 326, 327].

### System Specifications
The radar parameters were calculated to meet specific automotive constraints:

| Parameter | Value | Design Goal |
| :--- | :--- | :--- |
| **Carrier Frequency ($f_c$)** | 76.5 GHz | [cite_start]Automotive Standard [cite: 333] |
| **Bandwidth ($BW$)** | 1 GHz | [cite_start]Range Res $\le 0.75m$ [cite: 334] |
| **Range Resolution** | **0.15 m** | [cite_start]High Precision [cite: 336] |
| **Max Range ($R_{max}$)** | **315 m** | [cite_start]$> 250m$ Requirement [cite: 335] |
| **Max Velocity ($V_{max}$)** | **466.85 m/s** | [cite_start]$> 100m/s$ Requirement [cite: 335] |

---

## Signal Processing Chain

### 1. Waveform Generation
[cite_start]The system generates LFM chirps with a sweep time of $2.1 \mu s$[cite: 335].
![Transmitted Signal](images/TX.png)
*Figure 1: The transmitted FMCW sawtooth waveform and pulse train.*

### 2. Received Signal & Mixing
Echoes are received with time delays (range) and phase shifts (Doppler). [cite_start]White Gaussian Noise is added to simulate realistic SNR conditions (5dB, 10dB, 15dB)[cite: 446].
![Received Signal](images/RX.png)
*Figure 2: The received radar echo signal with additive noise.*

### 3. Detection & Range-Doppler Mapping
The beat signal is processed using a 2D-FFT. The 3D Range-Doppler Map visualizes targets as distinct peaks in the Range-Velocity plane.
![Range Doppler Map](images/range_doppler_map.png)
*Figure 3: 3D Range-Doppler Map showing clear target separation.*

---

## Performance Results

The algorithm was validated against "Ground Truth" targets. The detection error was negligible, demonstrating the robustness of the peak detection algorithm.

![Detection Results](images/detection_results.png)
*Figure 4: MATLAB Console output showing precise estimation of range and velocity.*

**Accuracy Summary:**
* [cite_start]**Range Error:** $< 0.03$ meters [cite: 823, 826]
* [cite_start]**Velocity Error:** $< 0.04$ m/s [cite: 833, 837]
* [cite_start]**Direction:** Successfully classified targets as "Approaching" or "Receding" [cite: 816-818].


## How to Run
1.  Clone this repository.
2.  Open `src/FMCW.m` in MATLAB.
3.  Run the script to generate the waveforms and detection results.
