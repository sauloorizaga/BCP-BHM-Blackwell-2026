# Disturbing the Peace at HPC Centers: Bringing 40x Data-Center Performance to the Individual Workstation.

# BCP-BHM-Blackwell-2026
First benchmark of 512³ 3D phase-field simulations  on consumer NVIDIA Blackwell GPUs via the BHM scheme  in MATLAB 2026a.  


# Why BHM?
Unlike SAV or Convex-Splitting methods, the Biharmonic-Modified (BHM) scheme avoids system enlargement and auxiliary variables. This linear, single-solve structure is essential for memory-dominated regimes, enabling 
 double-precision simulations on 24GB GPUs where other methods would typically trigger VRAM overflow.

Companion code for:
Orizaga et al. (2026), *Computers and Mathematics 
with Applications* (under review).

## Requirements
- MATLAB 2026a with Parallel Computing Toolbox
- NVIDIA GPU (RTX 5090 recommended)
- 24GB VRAM for N=512³ simulations

## Contents
- `CH3D_GPU_project_2025_BCP.m`: Main BCP simulation
- `visualization.m` — Isosurface and cross-section plots

## Usage
```matlab
% Run simulation
CH3D_GPU_project_2025_BCP(dt, M1, iter, tfinal, N)

% Example parameters
dt=0.01; M1=7.5; tfinal=100; N=256;
```
## Output
Running the codes: **CH3D_GPU_project_2025_BCP.m** followed by **visualization.m** will generate the following plot

<img src="BCP.png" width=400px height=400px>


## Citation
If you use this code please cite:
Orizaga et al. (2026), Computers and Mathematics 
with Applications.
