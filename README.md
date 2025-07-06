# Desktop-Software zur Ästuarzirkulation-Modellierung (Software for Estuarine Circulation Modelling)

## Funktionalitäten (Functionalities)

The Estuarine Circulation Modeling software is a desktop software written to simulate and analyze estuarine hydrodynamics, stratification, turbulence, sediment transport, tidal dynamics, and wave-current interactions. It integrates numerical solvers, turbulence models, and visualization tools to model complex estuarine processes. The software supports both structured and unstructured grids, multiple coordinate systems (sigma and z-level), and advanced numerical schemes like Total Variation Diminishing (TVD).

1. **Core Estuarine Model**:
   - Manages core estuarine parameters: estuary length (1000 m), depth (10 m), tidal amplitude (1 m), tidal period (12 hours), salinity (river: 0 PSU, ocean: 35 PSU), and temperature (river: 25°C, ocean: 20°C).
   - Supports hydrostatic and non-hydrostatic modes with an optional Reynolds-Averaged Navier-Stokes (RANS) solver.
   - Tracks salt wedge position and provides methods to retrieve salinity, temperature, and velocity profiles at specific points.

2. **Hydrodynamic Solver**:
   - Solves 2D shallow water dynamics for water level and velocities on a 100x100 grid.
   - Incorporates tidal forcing, wind stress, Coriolis force, and bottom friction.
   - Applies boundary conditions for river (x=0) and ocean (x=1000 m) interfaces.
   - Integrates wet/dry dynamics and wave-current interactions.

3. **2D Shallow Water Equations**:
   - Implements a 2D shallow water model with a GUI for visualizing water level, velocity, and salinity fields.
   - Allows user control over tidal amplitude, period, wind speed, wind direction, wave height, and grid size.
   - Supports wet/dry transitions and wave-enhanced friction via external modules.

4. **Total Variation Diminishing**:
   - Uses a TVD scheme with Harten-Lax-van Leer (HLL) flux to advect salinity, temperature, turbulent kinetic energy (k), and dissipation rate (epsilon) or specific dissipation rate (omega).
   - Supports structured (200x100) and unstructured (triangular mesh, 40x25 nodes) grids with bathymetry (shallower near river, x < 200 m).
   - GUI offers visualization modes: plan view, cross-section, contour, and quiver plots, with controls for depth, x-position, tidal period, and turbulence model (k-epsilon or k-omega).
   - Tracks fields: salinity, temperature, velocity (x, z), density, k, epsilon/omega, eddy viscosity, and diffusivity.

5. **Stratification**:
   - Models 1D vertical stratification of density, salinity, temperature, and passive scalars on a 100-point grid.
   - Computes gradient Richardson number and adjusts eddy viscosity using k-epsilon, k-omega, or constant turbulence models.
   - GUI allows control of mixing coefficient, critical Richardson number, and river scalar concentration, with visualization of density and scalar profiles.

6. **Baroclinic Flow**:
   - Simulates buoyancy-driven flows due to density gradients from salinity and temperature variations.
   - Updates velocity fields using finite differences, incorporating buoyancy effects.

7. **Passive Scalar Transport Equation**:
   - Solves 1D advection-diffusion for passive scalars (e.g., pollutants).
   - Integrates with baroclinic flow for consistent transport dynamics.

8. **Comprehanesive Forcing Mechanism**:
   - Combines tidal, wind, and wave forcing for hydrodynamic simulations.
   - GUI allows adjustment of forcing parameters (tidal amplitude, wind speed, wave height) and visualization of velocity and water level fields.

9. **Wave Current Interaction**:
   - Computes Stokes drift velocities (u, v) using linear wave theory for shallow water.
   - Calculates wave-enhanced bottom friction via a simplified Grant-Madsen model.
   - Supports dynamic updates of wave height, direction, period, and depth.

10. **Wind Forcing**:
    - Computes wind drag coefficient and stress components (tauX, tauY) based on wind speed, direction, and wave height.
    - Uses a parameterized drag coefficient dependent on wind speed and wave height.

11. **Wet & Dry Algorithm**:
    - Manages wet/dry transitions for tidal flats by tracking cell status based on a minimum depth threshold (Dmin).
    - Sets velocities, water level, and salinity to zero in dry cells and adjusts fluxes at wet/dry interfaces.
    - Ensures mass conservation through flux divergence corrections.

12. **Simpson-Hunter Mechanism Parameterization**:
    - Models tidal straining and internal tide effects in a sigma-coordinate system.
    - Computes Stokes drift, vertical velocity, and turbulent kinetic energy (TKE) production due to tidal shear.

13. **Asymmetrical Tidal Mixing**:
    - Simulates asymmetric tidal mixing effects on salinity and velocity fields.
    - Uses a `Cell` class to manage local properties (salinity, velocity, TKE).

14. **Bifurcated Estuary Model**:
    - Models circulation in bifurcated estuarine channels, accounting for branching flow dynamics and tidal effects.

15. **Equations of State**:
    - Computes water density based on salinity and temperature using equations of state.

16. **Vertical Discretization**:
    - Provides sigma (terrain-following) and z-level (fixed horizontal layers) coordinate systems for vertical grid setup.
    - Computes metric terms (z_xi, z_eta) for grid transformations.

17. **Richardson Number Dependent And SSI Mixing**:
    - Computes gradient Richardson number and strain-induced mixing for turbulence closure.
    - Supports k-epsilon and k-omega models for eddy viscosity calculations.

18. **Large Eddies Simulationcs**:
    - Implements Large Eddy Simulation (LES) for high-resolution turbulence modeling.
    - UI allows control of Smagorinsky coefficient and visualization of velocity and vorticity fields.

19. **Lattice-Boltzmann LES**:
    - Uses Lattice Boltzmann Method (LBM) for LES, focusing on vorticity and velocity fields.
    - UI supports adjustment of grid size and relaxation time, with visualization of vorticity contours.

20. **Spectral Analyzer**:
    - Performs spectral analysis (Welch’s method) and Empirical Orthogonal Function (EOF)/Proper Orthogonal Decomposition (POD) on variables like Richardson number, velocity, and salinity.
    - UI offers controls for variable selection, window type (Hanning, Hamming, Blackman, Rectangular), and visualization (PSD, heatmaps, spatial/temporal modes).

21. **Multi-Fraction Sediment Transporation**:
    - Models multi-fraction sediment transport, including bedload and suspended load.
    - UI allows adjustment of sediment properties (grain size, settling velocity) and visualization of concentration profiles.

22. **Adaptive Mesh Refinement**:
    - Implements adaptive mesh refinement based on velocity gradients or water level changes.
    - Dynamically adjusts grid resolution to capture fine-scale features.

23. **Visualization Renderer**:
    - Renders 2D visualizations of water level, salt wedge, salinity, temperature, passive scalar, and velocity profiles.
    - Used by UI modules for consistent plotting across simulations.

## Arbeitslogik (Operation Logic)

1. **Grid and Field Initialization**:
   - Modules like `EstuarineModel.cs`, `VerticalDiscretization.cs`, and `TotalVariationDiminishing.cs` initialize grids:
     - Structured grids (e.g., 200x100 in `TotalVariationDiminishing.cs`, 100x100 in `HydrodynamicSolver.cs`).
     - Unstructured triangular meshes (e.g., 40x25 nodes in `TotalVariationDiminishing.cs`).
     - Sigma or z-level coordinates (`VerticalDiscretization.cs`) account for bathymetry (shallower near river, x < 200 m).
   - Fields (salinity, temperature, velocity, density, TKE, epsilon/omega) are initialized with realistic profiles, e.g., salinity varying smoothly from 0 PSU at the river to 35 PSU at the ocean.

2. **Forcing and Boundary Conditions**:
   - Tidal forcing is applied via sinusoidal water level variations in `EstuarineModel.cs`, `ShallowWaterEq2D.cs`, and `TotalVariationDiminishing.cs`.
   - Wind forcing (`WindForcing.cs`) and wave-current interactions (`WaveCurrentInteraction.cs`) add surface and bottom stresses.
   - River (x=0, low salinity/temperature) and ocean (x=1000 m, high salinity/temperature) boundary conditions are set in `TotalVariationDiminishing.cs` and `Stratification.cs`.

3. **Numerical Solvers**:
   - `HydrodynamicSolver.cs` and `ShallowWaterEq2D.cs` use semi-implicit finite differences to solve shallow water dynamics for water level and velocities, with small time steps (e.g., 0.1 s) for stability.
   - `TotalVariationDiminishing.cs` uses a TVD scheme with HLL flux and a limiter to ensure stable advection of scalars and turbulence variables.
   - `LargeEddySim.cs` uses a Smagorinsky subgrid model for LES, while `LatticeBoltzmannLES.cs` employs LBM with a D3Q19 lattice for turbulence.
   - `PassiveScalarTransportEq.cs` and `Stratification.cs` solve advection-diffusion for scalars using finite differences.

4. **Turbulence Modeling**:
   - k-epsilon and k-omega models in `TotalVariationDiminishing.cs`, `Stratification.cs`, and `RichardsonNumDepAndSSIMix.cs` compute turbulent kinetic energy and dissipation/specific dissipation.
   - Eddy viscosity and diffusivity are dynamically calculated, modulated by shear and buoyancy effects via the gradient Richardson number.

5. **Transport and Stratification**:
   - `PassiveScalarTransportEq.cs` and `Stratification.cs` advect and diffuse scalars (e.g., pollutants) using velocity fields from `BaroclinicFlow.cs` or `HydrodynamicSolver.cs`.
   - Density gradients drive baroclinic flows in `BaroclinicFlow.cs`, coupled with density calculations from `EqOfState.cs`.

6. **Wet/Dry Dynamics**:
   - `WetAndDryAlgo.cs` updates cell status (wet/dry) based on total depth and adjusts fluxes to prevent flow into dry cells.
   - Ensures mass conservation by correcting water level using flux divergence.

7. **Visualization and Analysis**:
   - GUI modules (`TotalVariationDiminishing.cs`, `Stratification.cs`, `ShallowWaterEq2D.cs`, etc.) use `VisualizationRenderer.cs` to plot:
     - Plan view: Salinity, temperature, velocity magnitude at a selected depth.
     - Cross-section: Vertical profiles at a selected x-position.
     - Contour: Field distributions (structured: rectangular cells, unstructured: triangular mesh).
     - Quiver: Velocity vectors.
   - `SpectralAnalyzer.cs` computes power spectral density (PSD) and EOF/POD modes for diagnostic analysis of variables like velocity and salinity.

8. **Dynamic Updates**:
   - `AdaptiveMeshRef.cs` refines grids based on velocity gradients or water level changes.
   - User inputs (e.g., tidal period, wind speed) are updated via GUI controls or programmatic methods (e.g., `UpdateParameters` in `WindForcing.cs`).

## Physikalische und Mathematische Modelle (Physics and Mathematical Models)

1. **Shallow Water Dynamics**
   - Models water level and velocity fields in 2D, accounting for tidal forcing, wind stress, Coriolis force, and bottom friction. Includes mass conservation and momentum balance, with support for wet/dry transitions and wave effects.

2. **Baroclinic Flow**
   - Simulates buoyancy-driven flows caused by density gradients due to salinity and temperature variations. Density is computed using equation of state from `EqOfState.cs`, influencing velocity updates.

3. **Passive Scalar Transport**
   - Models the advection and diffusion of passive scalars (e.g., pollutants) in 1D, driven by velocity fields and eddy diffusivity. Includes source terms for river and ocean inputs.

4. **Turbulence Models**
   - **k-epsilon Model**: Computes turbulent kinetic energy (TKE) and dissipation rate, with shear and buoyancy production terms. Eddy viscosity is derived from TKE and dissipation.
   - **k-omega Model**: Computes TKE and specific dissipation rate as an alternative, with similar production terms. Eddy viscosity is based on TKE and specific dissipation.
   - **Gradient Richardson Number**: Modulates mixing based on stratification and shear, used to adjust eddy viscosity and diffusivity.

5. **Wave-Current Interaction**
   - Calculates Stokes drift velocities using linear wave theory for shallow water, accounting for wave height, period, and depth. Computes wave-enhanced bottom friction using a simplified Grant-Madsen model, adjusting the friction coefficient based on wave orbital velocity.

6. **Wind Forcing**
   - Computes wind stress components based on wind speed, direction, and wave height. Uses a parameterized drag coefficient that varies with wind speed and sea state, bounded for physical realism.

7. **Wet/Dry Algorithm**
   - Determines wet or dry status of grid cells based on a minimum depth threshold, accounting for water level and bathymetry. Adjusts velocities and scalars in dry cells and ensures mass conservation through flux corrections.

8. **Spectral Analysis**
   - Performs spectral analysis using Welch’s method to compute power spectral density for variables like velocity and salinity. Applies EOF/POD analysis via singular value decomposition to extract spatial and temporal modes.

9. **Sediment Transport**
   - Models multi-fraction sediment transport, including bedload (driven by bed shear stress) and suspended load (advected and diffused with settling velocity). Supports multiple grain sizes and dynamic sediment properties.

10. **Tidal Straining**
    - Models tidal straining and internal tide effects, computing Stokes drift and TKE production due to tidal shear in a sigma-coordinate system.
      
## Struktur der Ästuarzirkulation-Modelle (Estuarine Circulation Modeling Framework)

### Hydrodynamic Solver

The momentum equation solves:

```
∂u/∂t + u·∂u/∂x = -1/ρ ∂p/∂x + ν_eff ∂²u/∂x² + g ∂η/∂x + F_coriolis
```

Where:
- u = velocity (m/s)
- p = pressure (Pa)
- ρ = density (kg/m³)
- ν_eff = effective viscosity (m²/s)
- η = water surface elevation (m)
- F_coriolis = Coriolis force term

For non-hydrostatic cases, the pressure splits into:
p = p_hydrostatic + p_dynamic
where p_hydrostatic = ρ*g*(D - z)

### k-ε Turbulence Modeling

The RANS solver implements:
1. Turbulent kinetic energy (k) equation:
```
∂k/∂t + u·∂k/∂x = ∂/∂x[(ν+ν_t/σ_k)∂k/∂x] + P_k - ε
```
2. Dissipation rate (ε) equation:  
```
∂ε/∂t + u·∂ε/∂x = ∂/∂x[(ν+ν_t/σ_ε)∂ε/∂x] + C1_ε(ε/k)P_k - C2_ε(ε²/k)
```
With:
- ν_t = Cμ k²/ε (eddy viscosity)
- P_k = ν_t (∂u/∂x)² (production term)
- Constants: Cμ=0.09, σ_k=1.0, σ_ε=1.3, C1_ε=1.44, C2_ε=1.92

### Salinity/Temperature Transport

Solved via advection-diffusion:
```
∂C/∂t + u·∂C/∂x = ∂/∂x(K ∂C/∂x)
```
Where:
- C = scalar concentration (salinity or temperature)
- K = diffusivity (m²/s) = molecular + turbulent (ν_t/Sc_t)
- Sc_t = turbulent Schmidt number (~0.7)

### Bifurcation Handling

The bifurcated channel model:
1. Creates three connected segments:
   - Main channel (0 ≤ x ≤ L_bifurcation)
   - Branch 1 (scaled depth × 0.8)
   - Branch 2 (scaled depth × 1.2)
2. Modifies mixing at branches via:
   ```
   ν_t_branch = ν_t_main × depth_scale
   ```
   ```
   u_tidal_branch = u_tidal_main × tidal_scale
   ```
3. Maintains mass continuity at junction

### Numerical Model

- Spatial discretization: Central differences for diffusion terms
- Time integration: Explicit Euler scheme
- Stability enforcement:
  - Velocity CFL condition: u·dt/dx < 1
  - Diffusion condition: K·dt/dx² < 0.5
- Boundary conditions:
  - River end: Fixed inflow velocity
  - Ocean end: Tidal forcing
  - Branch ends: Zero-gradient

### Shallow Water Equations (Hydrodynamics)

**Continuity Equation**:
```
∂η/∂t + ∂(H·u)/∂x = 0
```

**Momentum Equation**:
```
∂u/∂t + u·∂u/∂x = -g ∂η/∂x + (1/H) ∂(H·τ)/∂x - (C_f·|u|·u)/H + F
```

Where:
- η = water surface elevation (m)
- H = total water depth = h_bottom + η
- u = depth-averaged velocity (m/s)
- τ = turbulent stress tensor (τ = ν_eff·∂u/∂x)
- C_f = bottom friction coefficient (~0.002-0.005)
- F = body forces (Coriolis, wind, baroclinic)

### Turbulence Closure (k-ε Model)

**Turbulent Kinetic Energy (k)**:
```
∂k/∂t + u·∂k/∂x = P_k + ∂/∂x[(ν+ν_t/σ_k)∂k/∂x] - ε
```

**Dissipation Rate (ε)**:
```
∂ε/∂t + u·∂ε/∂x = C1_ε(ε/k)P_k + ∂/∂x[(ν+ν_t/σ_ε)∂ε/∂x] - C2_ε(ε²/k)
```
- With production term: `P_k = ν_t·(∂u/∂x)²`
- Eddy viscosity calculated as: `ν_t = Cμ·k²/ε`

### Salinity Transport Equation
```
∂S/∂t + u·∂S/∂x = ∂/∂x[(K_m + ν_t/Sc_t)∂S/∂x]
```
Where:
- S = salinity (PSU)
- K_m = molecular diffusivity (~10⁻⁹ m²/s)
- Sc_t = turbulent Schmidt number (~0.7)

### Temperature Transport Equation
```
∂T/∂t + u·∂T/∂x = ∂/∂x[(α + ν_t/Pr_t)∂T/∂x] + Q_net/(ρ·c_p·H)
```
Where:
- T = temperature (°C)
- α = thermal diffusivity (~1.4×10⁻⁷ m²/s)
- Pr_t = turbulent Prandtl number (~0.85)
- Q_net = net heat flux (W/m²)
- c_p = heat capacity (~4186 J/kg·K)

### Salt Wedge Dynamics

The salt wedge position (x_w) evolves via: 
```
dx_w/dt = u_r - u_e - u_t·sin(ωt)
```
- u_r = river velocity = Q_r/(B·H)
- u_e = entrainment velocity = E·√(g'·H)
- u_t = tidal velocity amplitude
- g' = reduced gravity = g·Δρ/ρ
- E = entrainment coefficient (~10⁻⁴-10⁻³)

### Baroclinic Pressure Gradient

Implemented as:
```
∂p_b/∂x = -g ∫_η^-h (∂ρ/∂x) dz
```

With density calculated via:
```
ρ = ρ_fresh + β_S·S + β_T·T
```
where β_S ≈ 0.8 kg/m³/PSU, β_T ≈ -0.2 kg/m³/°C

### Bifurcation Physics

**Branch Flow Partitioning** uses:
```
Q_1/Q_2 = (B_1·H_1^(5/3)·√(S_1)) / (B_2·H_2^(5/3)·√(S_2))
```
- Q = discharge
- B = width
- H = depth
- S = energy slope

**Branch-Specific Modifications**:
```
ν_t^branch = ν_t^main × (H_branch/H_main)^(4/3)
```
```
u_tidal^branch = u_tidal^main × (H_main/H_branch)^(1/2)
```

### Numerical Discretization

**Spatial Derivatives**:
- Advection: 
```
∂φ/∂x ≈ (φ_{i+1} - φ_{i-1})/(2Δx)
```
- Diffusion: 
```
∂²φ/∂x² ≈ (φ_{i+1} - 2φ_i + φ_{i-1})/Δx²
```
**Time Integration**:
```
φ^{n+1} = φ^n + Δt·F(φ^n)
```
**Stability Criteria**:
- CFL: Δt < Δx/max(|u| + √(gH))
- Diffusion: Δt < 0.5·Δx²/ν_eff

### Boundary Conditions
1. **River Boundary (x=0)**:
   - u = Q_r/(B·H)
   - S = 0 PSU
   - T = T_river
2. **Ocean Boundary (x=L)**:
   - η = A_tide·sin(ωt)
   - S = S_ocean
   - T = T_ocean
3. **Branch Connections**:
   - Continuity: ΣQ_in = ΣQ_out
   - Momentum: p_1 = p_2 = p_3 at junction

## Gesamtvariation Abnehmend (Total Variation Diminishing)

### Fluid Dynamics Core

**Shallow Water Momentum Equations**:
```
∂u/∂t + u·∂u/∂x + w·∂u/∂z = -1/ρ ∂p/∂x + ∂/∂z(ν_t ∂u/∂z) + g_x
```
```
∂w/∂t + u·∂w/∂x + w·∂w/∂z = -1/ρ ∂p/∂z + ∂/∂z(ν_t ∂w/∂z) + g_z
```
- u = horizontal velocity (m/s)
- w = vertical velocity (m/s)
- ν_t = turbulent eddy viscosity (m²/s)
- g_x, g_z = gravitational acceleration components

**Density Calculation**:
```
ρ = ρ_0 [1 - α(T-T_0) + β_S(S-S_0)]
```
- α = thermal expansion coeff (2e-4 /°C)
- β_S = saline contraction coeff (8e-4 /PSU)

### Scalar Transport

**Salinity/Temperature Advection**:
```
∂C/∂t + ∇·(uC) = ∇·(K∇C) + S_C
```
Implemented via TVD flux limiter:
```
F_i+1/2 = F_L + ψ(r)(F_H - F_L)
```
- `ψ(r) = max(0, min(2r, (1+r)/2, 2))` (Superbee limiter)
- `r = (C_i - C_i-1)/(C_i+1 - C_i)` (gradient ratio)

### Turbulence Closure Models

**k-epsilon Equations**:
```
∂k/∂t + u·∇k = ∇·[(ν+ν_t/σ_k)∇k] + P_k - ε
```
```
∂ε/∂t + u·∇ε = ∇·[(ν+ν_t/σ_ε)∇ε] + C1_ε(ε/k)P_k - C2_ε(ε²/k)
```
**k-omega Equations**:
```
∂k/∂t + u·∇k = ∇·[(ν+ν_t/σ_k)∇k] + P_k - β* kω
```
```
∂ω/∂t + u·∇ω = ∇·[(ν+ν_t/σ_ω)∇ω] + α_ω(ω/k)P_k - β ω²
```
- `P_k = ν_t(∂u/∂z)² - (g/ρ_0)ν_t/Pr_t ∂ρ/∂z` (shear+buoyancy production)
- `ν_t = Cμ k²/ε` (k-epsilon) or `k/ω` (k-omega)

### TVD Flux Calculation

The HLL (Harten-Lax-van Leer) approximate Riemann solver:
```
F_HLL = 
{ F_L if s_L ≥ 0
{ F_R if s_R ≤ 0
{ (s_R F_L - s_L F_R + s_L s_R (U_R - U_L))/(s_R - s_L) otherwise
```
- `s_L = min(u_L - c_L, u_R - c_R, 0)`
- `s_R = max(u_L + c_L, u_R + c_R, 0)`

### Unstructured Grid Handling

For triangular elements:
- Gradient reconstruction via Green-Gauss method:
  ```
  ∇φ ≈ 1/V ∑_faces (φ_face · n_face) A_face
  ```
- Flux integration using control volume approach

### Stability Enforcement

- CFL condition: 
   ```
   max(|u|Δt/Δx + |w|Δt/Δz) ≤ 1
   ```
- Positivity preservation for k, ε, ω via:
  ```
  φ_new = max(φ_min, φ + Δt·RHS)
  ```
## Stratifizierung (Stratification)

**Key Components**:
- Salinity and temperature-driven density variations
- Turbulence closure schemes (k-epsilon/k-omega)
- Passive scalar transport
- Richardson number-based mixing suppression

### Density Calculation
Uses equation of state:
```
ρ = ρ₀[1 - α(T-T₀) + βₛ(S-S₀)]
```
- ρ = water density (kg/m³)
- α = thermal expansion coeff (2×10⁻⁴/°C)
- βₛ = saline contraction coeff (8×10⁻⁴/PSU)
- T₀,S₀ = reference values (20°C, 35PSU)

### Gradient Richardson Number
```
Ri = (g/ρ₀)(∂ρ/∂z)/(∂u/∂z)²
```
Critical values:
- `Ri > 0.25` → Strong stratification (reduced mixing)
- `Ri < 0.25` → Well-mixed conditions

### Turbulence Modulation

**k-epsilon Model**:
```
ν_eff = ν₀/(1 + Ri/Ri_crit)
```
**k-omega Model**:
```
ν_eff = ν₀(1 + 0.5 min(Ri,Ri_crit))/(1 + Ri)
```
- ν₀ is the neutral eddy viscosity

### Scalar Transport Equation
```
∂C/∂t + u·∂C/∂x = ∂/∂x(K ∂C/∂x) + S
```
With boundary conditions:
- River end: C = C_river
- Ocean end: C = C_ocean

### Mechanism
1. Compute density profile from S,T fields
2. Calculate Richardson number at each node
3. Adjust eddy viscosity based on Ri
4. Solve transport equations:
   - Salinity: `∂S/∂t + u·∇S = ∇·(K∇S)
   - Temperature: ∂T/∂t + u·∇T = ∇·(K∇T)
   - Passive scalar: ∂C/∂t + u·∇C = ∇·(K∇C)

### Stability Controls
- Density clamping: `1000 ≤ ρ ≤ 1030 kg/m³`
- Scalar bounds: `0 ≤ C ≤ C_river`
- Time step limitation: `Δt < Δx²/(2K_max)`

## Baroklinische-Strömung (Baroclinic Flow)

**Key Components**:
- Equation of state for seawater density
- Richardson number-dependent turbulent mixing
- Heat flux parameterization (shortwave, longwave, sensible, latent)
- Semi-implicit numerical schemes for stability

### Density Calculation
Equation of state:
```
ρ(S,T,P) = ρ₀ + Δρ(S,T,P)
```
- ρ₀ = 1000 kg/m³ (reference)
- Δρ = density anomaly from salinity (S), temperature (T), pressure (P)

### Baroclinic Pressure Gradient
```
∇p_b = g·h·∇ρ/ρ₀
```
- g = 9.81 m/s² (gravity)
- h = estuary depth (m)
- ∇ρ = horizontal density gradient

### Turbulent Mixing

**Eddy Diffusivities**:
- Horizontal: 
  ```
  Kₓ = νₜ + l·√(|∂u/∂x|)
  ```
- Vertical: 
  ```
  K_z = (νₜ + K_z₀)/(1 + Ri/Ri_c)
  ```
   - νₜ = turbulent viscosity (m²/s)
   - l = 10m (mixing length scale)
   - Ri = g/ρ₀·∂ρ/∂z/(∂u/∂z)² (Richardson number)
   - Ri_c = 0.25 (critical value)

### Heat Flux Components

**Surface Energy Balance**:
```
Q_net = Q_sw↓ - Q_lw↑ - Q_sensible - Q_latent
```
- Shortwave: `Q_sw = 1000(1-α)(1-cc)cos(2πt/86400)` W/m²
- Longwave: `Q_lw = εσ(T_w⁴ - T_a⁴)`
- Sensible: `Q_sen = ρₐc_pₐC_hU(T_w-T_a)B`
- Latent: `Q_lat = ρₐL_vC_eU(q_s-q_a)`

### Transport Equations

**Salinity**:
```
∂S/∂t + u·∇S = ∇·(Kₓ∂S/∂x + K_z∂S/∂z) + Q_S
```
- Q_S = -γ·|∇u|·S (river) or +γ·|∇u|·(S_ocean-S) (ocean)

**Temperature**:
```
∂T/∂t + u·∇T = ∇·(Kₓ∂T/∂x + K_z∂T/∂z) + Q_T/ρc_p
```
- Q_T = net surface heat flux

### Mechanism
1. **Advection**:
   - Van Leer flux limiter scheme
   - 2nd order accurate with TVD properties
2. **Diffusion**:
   - Crank-Nicolson semi-implicit method
   - Tridiagonal matrix solver (Thomas algorithm)
3. **Boundary Conditions**:
   - River: S=0, T=15°C
   - Ocean: S=35PSU, T=20°C

### Key Parameters

**Physical Constants**:
- Gravitational acceleration: 9.81 m/s²
- Water properties:
  - ρ₀ = 1000 kg/m³
  - c_p = 4184 J/kg·K
  - α = 2×10⁻⁴ /°C
  - β_S = 8×10⁻⁴ /PSU
- **Atmospheric Forcing**:
   - Wind speed: 5 m/s
   - Air temperature: 15°C
   - Cloud cover: 50%
   - Albedo: 6%
- **Numerical Settings**:
   - Base diffusivities:
   - Kₓ₀ = 0.1 m²/s
   - K_z₀ = 0.01 m²/s
   - Time step: 3600s (adaptive)
   - Grid resolution: Δx = L/N

### Implementation

1. **Stability Controls**:
   - Density clamping (1000-1030 kg/m³)
   - Temperature bounds (0-30°C)
   - Salinity bounds (0-35PSU)
2. **Performance Optimizations**:
   - Precomputed coefficients
   - Vectorized operations
   - Parallel tridiagonal solver
3. **Validation Metrics**:
   - Richardson number distribution
   - Salt wedge penetration
   - Surface heat budget closure

## Umfassender-Leistungmechanismus (Comprehensive Forcing Mechanism)

- **3D Navier-Stokes Solver**  
Solves the incompressible Navier-Stokes equations with Boussinesq approximation and salinity/temperature transport:

   - **Momentum Equations**:  
  ```
  ∂u/∂t + u·∇u = -∇p/ρ₀ + ν∇²u + g(ρ/ρ₀) + F_wind + F_Coriolis
  ```  
  (where ρ = ρ₀ + β_S·S - β_T·T)

   - **Continuity**: `∇·u = 0`  
   - **Salinity Transport**: `∂S/∂t + u·∇S = κ∇²S`  
   - **Pressure Correction**: Iteratively solves ∇²p = ∇·(u·∇u) to enforce incompressibility.

- **2D Shallow Water Equations**  
Depth-averaged form with wet/dry tracking:

   - **Mass Conservation**: `∂η/∂t + ∇·(H·u) = 0`  
  (H = h + η, where h is bathymetry, η is surface elevation)

   - **Momentum**:  
  ```
  ∂u/∂t + u·∇u = -g∇η + ν∇²u + τ_wind/(ρ₀H) - c_f|u|u/H  
  ```

- **1D Simplified Model**  
Longitudinal profiles only:
   - **Velocity**: 
   ```
   ∂u/∂t + u ∂u/∂x = -g ∂η/∂x + ν ∂²u/∂x² + τ_wind/(ρ₀D)  
   ```
   - **Salinity**: 
   ```
   ∂S/∂t + u ∂S/∂x = κ ∂²S/∂x²  
   ```
- **Flux Limiter Application** (Minmod limiter):  
   For a quantity φ (velocity, salinity):  
   ```
   φ_limited = φ + ψ(r)(φ_upwind - φ)  
   ```
   where ψ(r) = max(0, min(1, r)) and r = (φ - φ_upwind)/(φ_downwind - φ).
- **Advection Terms** (e.g., 3D salinity):  
   ```
   sAdvX = u ≥ 0 ? u·(S[i] - S[i-1])/Δx + ψ(r)(S[i+1] - S[i])/Δx  
           : u·(S[i+1] - S[i])/Δx + ψ(r)(S[i] - S[i-1])/Δx  
   ```
- **Time Step Control**:  
   ```
   Δt = CFL × min(Δx/(|u|+√gH), Δx²/(2ν))  
   ```
   (CFL = 0.4 for stability)
- **Tidal Forcing**:  
   ```
   η(t) = A_tide·sin(2πt/T_tide)  
   ```
   Implemented as boundary condition at the ocean end.
- **Wind Stress**:  
   ```
   τ_wind = ρ₀·c_d·|U_wind|²
   ```
   ```  
   c_d = 0.001, U_wind = (u_wind, v_wind)
   ```
- **River Discharge**:  
   Upstream boundary condition: u(0) = Q_river/(W·D).
- **Salinity Gradients**:  
   Source term in transport equation: ∇S = ∂S/∂x (user-defined).
- **Cell Status**:  
   Wet if H > h_min (default h_min = 0.01 m), else dry.
- **Boundary Handling**:  
   - No flux through dry cells.  
   - Momentum terms set to zero in dry regions.


