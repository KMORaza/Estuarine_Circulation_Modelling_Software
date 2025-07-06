using System;

namespace EstuarineCirculationModeling
{
    public class BaroclinicFlow
    {
        private readonly int gridPoints;
        private readonly double dx; // Spatial step
        private readonly double estuaryDepth;
        private readonly double densityReference = 1000.0; // Reference density (kg/m³)
        private readonly double g = 9.81; // Gravitational acceleration (m/s²)
        private readonly double baseHorizontalDiffusivity = 0.1; // Base Kx (m²/s)
        private readonly double baseVerticalDiffusivity = 0.01; // Base Kz (m²/s)
        private readonly double mixingRate = 0.005; // Mixing rate for source/sink terms (s^-1)
        private readonly double turbulentLengthScale = 10.0; // Turbulent length scale (m)
        private readonly double criticalRichardson = 0.25; // Critical Richardson number
        // Heat flux parameters
        private readonly double specificHeatCapacity = 4184.0; // Specific heat capacity of water (J/kg·K)
        private readonly double stefanBoltzmann = 5.67e-8; // Stefan-Boltzmann constant (W/m²·K⁴)
        private readonly double albedo = 0.06; // Surface albedo for water
        private readonly double emissivity = 0.97; // Emissivity of water
        private readonly double latentHeatVaporization = 2.45e6; // Latent heat of vaporization (J/kg)
        private readonly double bowenRatio = 0.61; // Bowen ratio for sensible heat (dimensionless)
        private readonly double windSpeed = 5.0; // Default wind speed (m/s)
        private readonly double cloudCover = 0.5; // Cloud cover fraction (0 to 1)
        private readonly double atmosphericTemperature = 15.0; // Atmospheric temperature (°C)

        public BaroclinicFlow(double estuaryLength, double estuaryDepth, int gridPoints)
        {
            this.gridPoints = gridPoints;
            this.dx = estuaryLength / gridPoints;
            this.estuaryDepth = estuaryDepth;
        }

        // Calculate density based on salinity and temperature
        private double CalculateDensity(double salinity, double temperature)
        {
            // Simplified equation of state: ρ = ρ₀ + 0.8S - 0.2T
            return densityReference + 0.8 * salinity - 0.2 * temperature;
        }

        // Compute baroclinic pressure gradient (∂p_b/∂x) at each grid point
        public double[] ComputeBaroclinicGradient(double[] salinity, double[] temperature)
        {
            double[] baroclinicGradient = new double[gridPoints];

            // Calculate density profile
            double[] density = new double[gridPoints];
            for (int i = 0; i < gridPoints; i++)
            {
                density[i] = CalculateDensity(salinity[i], temperature[i]);
            }

            // Compute baroclinic pressure gradient: ∂p_b/∂x ≈ g * ∫ (∂ρ/∂x) dz
            for (int i = 1; i < gridPoints - 1; i++)
            {
                // Approximate horizontal density gradient: ∂ρ/∂x
                double densityGradient = (density[i + 1] - density[i - 1]) / (2 * dx);
                // Baroclinic pressure gradient: g * ρ * h * (∂ρ/∂x) / ρ
                baroclinicGradient[i] = g * density[i] * estuaryDepth * densityGradient / densityReference;
            }

            // Boundary conditions
            baroclinicGradient[0] = baroclinicGradient[1];
            baroclinicGradient[gridPoints - 1] = baroclinicGradient[gridPoints - 2];

            return baroclinicGradient;
        }

        // Compute horizontal and vertical eddy diffusivities (Kx, Kz)
        private (double[] Kx, double[] Kz) ComputeEddyDiffusivities(double[] velocity, double[] salinity, double[] temperature, double[] eddyViscosity = null)
        {
            double[] Kx = new double[gridPoints];
            double[] Kz = new double[gridPoints];

            for (int i = 0; i < gridPoints; i++)
            {
                // Horizontal diffusivity: Kx = ν_t + l * sqrt(|∂u/∂x|)
                double nuT = eddyViscosity != null ? eddyViscosity[i] : baseHorizontalDiffusivity;
                double velocityGradient = i > 0 && i < gridPoints - 1
                    ? Math.Abs((velocity[i + 1] - velocity[i - 1]) / (2 * dx))
                    : 0.0;
                Kx[i] = nuT + turbulentLengthScale * Math.Sqrt(velocityGradient);

                // Vertical diffusivity: Kz = (ν_t + base) / (1 + Ri/Ri_c)
                double densityGradient = i > 0 && i < gridPoints - 1
                    ? (CalculateDensity(salinity[i + 1], temperature[i + 1]) - CalculateDensity(salinity[i - 1], temperature[i - 1])) / (2 * dx)
                    : 0.0;
                double shear = velocityGradient * velocityGradient;
                double Ri = shear > 1e-10 ? g * densityGradient / (densityReference * shear) : 0.0;
                Kz[i] = (nuT + baseVerticalDiffusivity) / (1 + Math.Max(0, Ri / criticalRichardson));
            }

            return (Kx, Kz);
        }

        // Compute surface heat flux (shortwave, longwave, sensible, latent)
        private double ComputeSurfaceHeatFlux(double waterTemp, double velocityGradient, double currentTime)
        {
            // Shortwave radiation: Q_sw = Q_solar * (1 - albedo) * (1 - cloudCover) * cos(ωt)
            double solarConstant = 1000.0; // Solar radiation at surface (W/m²)
            double diurnalCycle = Math.Cos(2 * Math.PI * currentTime / (24 * 3600)); // Diurnal cycle (24h period)
            double shortwaveFlux = solarConstant * (1 - albedo) * (1 - cloudCover) * Math.Max(0, diurnalCycle);

            // Longwave radiation: Q_lw = ε * σ * (T_w^4 - T_a^4)
            double waterTempK = waterTemp + 273.15; // Convert to Kelvin
            double atmTempK = atmosphericTemperature + 273.15;
            double longwaveFlux = emissivity * stefanBoltzmann * (Math.Pow(waterTempK, 4) - Math.Pow(atmTempK, 4));

            // Sensible heat: Q_sensible = ρ_air * C_p_air * C_h * U * (T_w - T_a)
            double airDensity = 1.225; // Air density (kg/m³)
            double airSpecificHeat = 1005.0; // Specific heat of air (J/kg·K)
            double sensibleHeatCoefficient = 0.0012; // Transfer coefficient
            double sensibleFlux = airDensity * airSpecificHeat * sensibleHeatCoefficient * windSpeed * (waterTemp - atmosphericTemperature);

            // Latent heat: Q_latent = ρ_air * L_v * C_e * U * (q_s - q_a)
            double latentHeatCoefficient = 0.0015; // Transfer coefficient
            double saturationVaporPressure = 6.1078 * Math.Pow(10, 7.5 * waterTemp / (waterTemp + 237.3)); // Tetens formula (hPa)
            double airVaporPressure = 6.1078 * Math.Pow(10, 7.5 * atmosphericTemperature / (atmosphericTemperature + 237.3)) * 0.8; // 80% humidity
            double latentFlux = airDensity * latentHeatVaporization * latentHeatCoefficient * windSpeed * (saturationVaporPressure - airVaporPressure) / 100.0; // Convert hPa to Pa

            // Total heat flux: Q_surface = Q_sw - Q_lw - Q_sensible - Q_latent (W/m²)
            double totalHeatFlux = shortwaveFlux - longwaveFlux - sensibleFlux - latentFlux;

            // Modulate by velocity gradient for mixing enhancement
            totalHeatFlux *= (1 + velocityGradient);

            // Convert to volumetric source term: Q_t = Q_surface / (ρ * Cp * h) (K/s)
            return totalHeatFlux / (densityReference * specificHeatCapacity * estuaryDepth);
        }

        // Solve salinity transport: ∂S/∂t + ∇·(uS) = ∇·(Kx·∂S/∂x + Kz·∂S/∂z) + Q_s
        public double[] SolveSalinityTransport(double[] salinity, double[] velocity, double dt, double salinityOcean, double saltWedgePosition, double estuaryLength, double[] eddyViscosity = null)
        {
            double[] newSalinity = new double[gridPoints];
            var (Kx, Kz) = ComputeEddyDiffusivities(velocity, salinity, new double[gridPoints], eddyViscosity);

            // Compute advective fluxes using Van Leer flux limiter
            double[] flux = new double[gridPoints + 1];
            for (int i = 0; i < gridPoints; i++)
            {
                double u = i < gridPoints - 1 ? (velocity[i] + velocity[Math.Min(i + 1, gridPoints - 1)]) / 2 : velocity[i];
                double sLeft = i > 0 ? salinity[i - 1] : salinity[0];
                double sRight = i < gridPoints - 1 ? salinity[i + 1] : salinity[gridPoints - 1];
                double sCurrent = salinity[i];

                // Van Leer limiter
                double r = u >= 0 ? (sCurrent - sLeft) / (sRight - sCurrent + 1e-10) : (sRight - sCurrent) / (sCurrent - sLeft + 1e-10);
                double phi = Math.Max(0, Math.Min(2 * r / (1 + r), 2));
                double sInterface = u >= 0 ? sCurrent + 0.5 * phi * (sRight - sCurrent) : sLeft + 0.5 * phi * (sCurrent - sLeft);

                flux[i] = u * sInterface;
            }
            flux[gridPoints] = velocity[gridPoints - 1] * salinity[gridPoints - 1];

            // Semi-implicit diffusion (Crank-Nicolson)
            double[] a = new double[gridPoints];
            double[] b = new double[gridPoints];
            double[] c = new double[gridPoints];
            double[] d = new double[gridPoints];

            for (int i = 1; i < gridPoints - 1; i++)
            {
                double x = i * dx;
                double kappaX = (Kx[i] + Kx[i + 1]) / 2; // Horizontal diffusion
                double kappaZ = Kz[i]; // Vertical diffusion
                double alphaX = kappaX * dt / (2 * dx * dx);
                double alphaZ = kappaZ * dt / (2 * estuaryDepth * estuaryDepth); // Parameterized vertical mixing

                // Source/sink term: River dilution and ocean mixing
                double velocityGradient = Math.Abs((velocity[Math.Min(i + 1, gridPoints - 1)] - velocity[Math.Max(i - 1, 0)]) / (2 * dx));
                double sourceSink = x < saltWedgePosition
                    ? -mixingRate * velocityGradient * salinity[i] // River dilution
                    : mixingRate * velocityGradient * (salinityOcean - salinity[i]); // Ocean mixing

                // Explicit advection and source terms
                double advection = (flux[i + 1] - flux[i]) / dx;
                double verticalMixing = alphaZ * (salinityOcean - salinity[i]); // Parameterized vertical diffusion
                d[i] = salinity[i] + dt * (-advection + sourceSink + verticalMixing);

                // Implicit horizontal diffusion
                a[i] = -alphaX;
                b[i] = 1 + 2 * alphaX;
                c[i] = -alphaX;
                d[i] += alphaX * (salinity[i + 1] - 2 * salinity[i] + salinity[i - 1]);
            }

            // Boundary conditions
            newSalinity[0] = 0.0; // River (freshwater)
            newSalinity[gridPoints - 1] = salinityOcean; // Ocean
            a[0] = 0; b[0] = 1; c[0] = 0; d[0] = newSalinity[0];
            a[gridPoints - 1] = 0; b[gridPoints - 1] = 1; c[gridPoints - 1] = 0; d[gridPoints - 1] = newSalinity[gridPoints - 1];

            // Solve tridiagonal system for diffusion
            newSalinity = SolveTridiagonal(a, b, c, d);

            // Ensure physical bounds
            for (int i = 0; i < gridPoints; i++)
            {
                newSalinity[i] = Math.Max(0, Math.Min(salinityOcean, newSalinity[i]));
            }

            return newSalinity;
        }

        // Solve temperature transport: ∂T/∂t + ∇·(uT) = ∇·(Kx·∂T/∂x + Kz·∂T/∂z) + Q_t
        public double[] SolveTemperatureTransport(double[] temperature, double[] velocity, double dt, double temperatureOcean, double saltWedgePosition, double estuaryLength, double[] eddyViscosity = null)
        {
            double[] newTemperature = new double[gridPoints];
            var (Kx, Kz) = ComputeEddyDiffusivities(velocity, new double[gridPoints], temperature, eddyViscosity);

            // Compute advective fluxes using Van Leer flux limiter
            double[] flux = new double[gridPoints + 1];
            for (int i = 0; i < gridPoints; i++)
            {
                double u = i < gridPoints - 1 ? (velocity[i] + velocity[Math.Min(i + 1, gridPoints - 1)]) / 2 : velocity[i];
                double tLeft = i > 0 ? temperature[i - 1] : temperature[0];
                double tRight = i < gridPoints - 1 ? temperature[i + 1] : temperature[gridPoints - 1];
                double tCurrent = temperature[i];

                // Van Leer limiter
                double r = u >= 0 ? (tCurrent - tLeft) / (tRight - tCurrent + 1e-10) : (tRight - tCurrent) / (tCurrent - tLeft + 1e-10);
                double phi = Math.Max(0, Math.Min(2 * r / (1 + r), 2));
                double tInterface = u >= 0 ? tCurrent + 0.5 * phi * (tRight - tCurrent) : tLeft + 0.5 * phi * (tCurrent - tLeft);

                flux[i] = u * tInterface;
            }
            flux[gridPoints] = velocity[gridPoints - 1] * temperature[gridPoints - 1];

            // Semi-implicit diffusion (Crank-Nicolson)
            double[] a = new double[gridPoints];
            double[] b = new double[gridPoints];
            double[] c = new double[gridPoints];
            double[] d = new double[gridPoints];

            for (int i = 1; i < gridPoints - 1; i++)
            {
                double x = i * dx;
                double kappaX = (Kx[i] + Kx[i + 1]) / 2; // Horizontal diffusion
                double kappaZ = Kz[i]; // Vertical diffusion
                double alphaX = kappaX * dt / (2 * dx * dx);
                double alphaZ = kappaZ * dt / (2 * estuaryDepth * estuaryDepth); // Parameterized vertical mixing

                // Surface heat flux
                double velocityGradient = Math.Abs((velocity[Math.Min(i + 1, gridPoints - 1)] - velocity[Math.Max(i - 1, 0)]) / (2 * dx));
                double heatFlux = ComputeSurfaceHeatFlux(temperature[i], velocityGradient, x / estuaryLength * 24 * 3600); // Scale time to diurnal cycle

                // Explicit advection and heat flux
                double advection = (flux[i + 1] - flux[i]) / dx;
                double verticalMixing = alphaZ * (temperatureOcean - temperature[i]); // Parameterized vertical diffusion
                d[i] = temperature[i] + dt * (-advection + heatFlux + verticalMixing);

                // Implicit horizontal diffusion
                a[i] = -alphaX;
                b[i] = 1 + 2 * alphaX;
                c[i] = -alphaX;
                d[i] += alphaX * (temperature[i + 1] - 2 * temperature[i] + temperature[i - 1]);
            }

            // Boundary conditions
            newTemperature[0] = 15.0; // River (fixed temperature)
            newTemperature[gridPoints - 1] = temperatureOcean; // Ocean
            a[0] = 0; b[0] = 1; c[0] = 0; d[0] = newTemperature[0];
            a[gridPoints - 1] = 0; b[gridPoints - 1] = 1; c[gridPoints - 1] = 0; d[gridPoints - 1] = newTemperature[gridPoints - 1];

            // Solve tridiagonal system for diffusion
            newTemperature = SolveTridiagonal(a, b, c, d);

            // Ensure physical bounds
            for (int i = 0; i < gridPoints; i++)
            {
                newTemperature[i] = Math.Max(0, Math.Min(temperatureOcean, newTemperature[i]));
            }

            return newTemperature;
        }

        // Solve tridiagonal system using Thomas algorithm
        private double[] SolveTridiagonal(double[] a, double[] b, double[] c, double[] d)
        {
            double[] result = new double[gridPoints];
            double[] cPrime = new double[gridPoints];
            double[] dPrime = new double[gridPoints];

            // Forward elimination
            cPrime[0] = c[0] / b[0];
            dPrime[0] = d[0] / b[0];

            for (int i = 1; i < gridPoints; i++)
            {
                double m = b[i] - a[i] * cPrime[i - 1];
                cPrime[i] = c[i] / m;
                dPrime[i] = (d[i] - a[i] * dPrime[i - 1]) / m;
            }

            // Back substitution
            result[gridPoints - 1] = dPrime[gridPoints - 1];
            for (int i = gridPoints - 2; i >= 0; i--)
            {
                result[i] = dPrime[i] - cPrime[i] * result[i + 1];
            }

            return result;
        }
    }
}