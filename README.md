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
   - Simulates buoyancy-driven flows caused by density gradients due to salinity and temperature variations. Density is computed using a linear or UNESCO equation of state from `EqOfState.cs`, influencing velocity updates.

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
