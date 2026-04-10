# BCP-BHM-Blackwell-2026
First benchmark of 512³ 3D phase-field simulations  on consumer NVIDIA Blackwell GPUs via the BHM scheme  in MATLAB 2026a. Companion code for Orizaga et al. (2026).

Companion code for:
Orizaga et al. (2026), *Computers and Mathematics 
with Applications* (under review).

## Requirements
- MATLAB 2026a with Parallel Computing Toolbox
- NVIDIA GPU (RTX 5090 recommended)
- 24GB VRAM for N=512³ simulations

## Contents
- `CH3D_GPU_project_2025_BCP.m` — Main BCP simulation
- `visualization.m` — Isosurface and cross-section plots

## Usage
```matlab
% Run simulation
CH3D_GPU_project_2025_BCP(dt, M1, iter, tfinal, N)

% Example parameters
dt=0.01; M1=7.5; tfinal=100; N=256;
```

## Citation
If you use this code please cite:
Orizaga et al. (2026), Computers and Mathematics 
with Applications.
