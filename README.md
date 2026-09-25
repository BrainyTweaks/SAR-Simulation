# SAR Simulation in MATLAB

This is a MATLAB SAR simulation project I made to understand the basic signal processing involved in Synthetic Aperture Radar imaging.

The simulation starts with an LFM chirp and generates echoes from a few simulated targets while the radar moves along a synthetic aperture. Noise is added to the received data, followed by range compression and SAR image formation using backprojection.

The main MATLAB file is `sar_simulation.m`.

## Radar Parameters

- Carrier frequency: $f_c = 5.3\,\text{GHz}$
- Bandwidth: $B = 40\,\text{MHz}$
- Pulse duration: $T_p = 5\,\mu\text{s}$
- Chirp rate: $K = \frac{B}{T_p} = 8.00 \times 10^{12}\,\text{Hz/s}$
- Platform velocity: $v = 100\,\text{m/s}$
- Center range: $R_0 = 5000\,\text{m}$
- Sampling frequency: $f_s = 80\,\text{MHz}$
- Synthetic aperture length: $L = 200\,\text{m}$
- Number of aperture positions: $N_a = 401$

The wavelength is:

$$
\lambda = \frac{c}{f_c} \approx 0.0566\,\text{m}
$$

The theoretical range resolution is:

$$
\Delta R = \frac{c}{2B} = 3.75\,\text{m}
$$

## Target Scene

Three point targets are used in the simulation.

| Target | X Position (m) | Range (m) | Reflectivity |
|---|---:|---:|---:|
| T1 | -40 | 5000 | 1.0 |
| T2 | 0 | 5000 | 0.7 |
| T3 | 40 | 5050 | 0.4 |

The different reflectivity values are used to check whether the relative target strengths are preserved after SAR processing.

## Transmitted LFM Chirp

The simulation starts by generating an LFM chirp.

The chirp rate is:

$$
K = \frac{B}{T_p}
$$

The transmitted signal is represented as:

$$
s(t) =
\exp\left[
j2\pi
\left(
f_c t + \frac{K}{2}t^2
\right)
\right]
$$

![Transmitted LFM Chirp](figures/01%20Transmitted%20LFM%20Chirp.png)

## Synthetic Aperture

The radar is simulated at different positions over a 200 m aperture.

![Synthetic Aperture](figures/02%20Synthetic%20Aperture.png)

## Raw SAR Data

For every radar position, the distance to each target is calculated and the corresponding delayed echo is generated. Noise is then added to the received signal.

![Raw SAR Data](figures/03%20Raw%20SAR%20Data.png)

## Range Compression

Range compression is performed using a matched filter for the transmitted chirp. FFT-based convolution is used for the filtering.

![Range-Compressed SAR Data](figures/04%20Range-Compressed%20SAR%20Data.png)

The theoretical range resolution for the 40 MHz bandwidth is:

$$
\Delta R = \frac{c}{2B} = 3.75\,\text{m}
$$

## SAR Focusing

The range-compressed data is focused using backprojection.

For each point in the image grid, the distance to every radar position is calculated. The corresponding samples from the range-compressed data are obtained using linear interpolation. The phase is then compensated and the signals from the aperture are coherently summed.

![Focused SAR Data](figures/05%20Focused%20SAR%20Data.png)

The resulting image contains the three simulated targets at approximately their expected positions.

## Focused SAR Image

The focused image is also displayed on a dB scale.

![Focused SAR Image](figures/06%20Focused%20SAR%20Image%20(dB).png)

The strongest response comes from T1, followed by T2 and T3, which matches their assigned reflectivities.

There are some sidelobe-like responses around the targets, but the three main targets are clearly localized.

## Range Profile

A range profile is taken through the strongest azimuth position.

![SAR Range Profile](figures/07%20SAR%20Range%20Profile.png)

## Azimuth Profile

An azimuth profile is extracted through the strongest range position.

![SAR Azimuth Profile](figures/08%20SAR%20Azimuth%20Profile.png)

## Peak Detection

The script also searches for the three strongest spatially separated image peaks.

The expected target locations are:

- T1: $X = -40\,\text{m}$, $R = 5000\,\text{m}$
- T2: $X = 0\,\text{m}$, $R = 5000\,\text{m}$
- T3: $X = 40\,\text{m}$, $R = 5050\,\text{m}$

The detected positions are printed in the MATLAB command window and can be compared with the known target positions.

## Running

The simulation was developed using MATLAB R2025b.

Open `sar_simulation.m` in MATLAB and run the script.

No external dataset is required. The target scene and received radar signals are generated directly by the simulation.

## Current Status

The basic SAR simulation is working, including:

- LFM chirp generation
- Target echo generation
- Noise addition
- Range compression
- Backprojection focusing
- Focused SAR image generation
- dB image display
- Range and azimuth profiles
- Basic peak detection

This is currently a simulation rather than a complete real-world SAR processing system. Further work can focus on sidelobe reduction, more realistic radar effects, and quantitative analysis of the focused image.