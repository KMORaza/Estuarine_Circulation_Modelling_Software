using System;
using System.Linq;

namespace EstuarineCirculationModeling
{
    public class ShallowWaterEq2D
    {
        public double[,] U, V, Eta, Salinity;
        private double estuaryLength, estuaryWidth, estuaryDepth;
        private int nx, ny;
        private double dx, dy;
        private double riverDischarge, tidalAmplitude, tidalPeriod, salinityGradient, waveHeight, stormSurgeAmplitude;
        private bool enableWetDry; // Added for wetting and drying
        private double minDepth;   // Added for minimum depth
        private readonly WindForcing windForcing; // Use WindForcing instance
        private double g = 9.81; // m/s²
        private double rho0 = 1000.0; // kg/m³ (freshwater density)
        private double rhoOcean = 1025.0; // kg/m³ (ocean water density)
        private double nu = 1e-6; // Kinematic viscosity (m²/s)
        private double kappa = 1e-4; // Salinity diffusion coefficient (m²/s)
        private const double courantNumber = 0.4; // CFL safety factor
        private const double frictionCoefficient = 0.0025; // Bottom friction coefficient
        private const double coriolisParameter = 1e-4; // s^-1 (mid-latitude approximation)
        private const double eddyViscosity = 0.01; // m²/s (turbulence closure)
        private const double atmosphericPressureGradient = 0.0001; // Pa/m
        private const double seasonalSalinityAmplitude = 2.0; // PSU

        public ShallowWaterEq2D(double length, double width, double depth, int nx, int ny,
            double riverDischarge, double tidalAmplitude, double tidalPeriod,
            WindForcing windForcing, double salinityGradient,
            double waveHeight, double stormSurgeAmplitude, bool enableWetDry, double minDepth)
        {
            this.estuaryLength = length;
            this.estuaryWidth = width;
            this.estuaryDepth = depth;
            this.nx = nx;
            this.ny = ny;
            this.riverDischarge = riverDischarge;
            this.tidalAmplitude = tidalAmplitude;
            this.tidalPeriod = tidalPeriod;
            this.windForcing = windForcing;
            this.salinityGradient = salinityGradient;
            this.waveHeight = waveHeight;
            this.stormSurgeAmplitude = stormSurgeAmplitude;
            this.enableWetDry = enableWetDry; // Store enableWetDry
            this.minDepth = minDepth;         // Store minDepth

            dx = estuaryLength / (nx - 1);
            dy = estuaryWidth / (ny - 1);
            U = new double[nx, ny];
            V = new double[nx, ny];
            Eta = new double[nx, ny];
            Salinity = new double[nx, ny];

            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    U[i, j] = riverDischarge / (estuaryWidth * estuaryDepth);
                    V[i, j] = 0.0;
                    Eta[i, j] = tidalAmplitude * Math.Sin(2 * Math.PI * (i * dx) / estuaryLength);
                    Salinity[i, j] = salinityGradient * i * dx;
                }
        }

        public double ComputeCFLTimeStep()
        {
            double maxU = U.Cast<double>().Max(Math.Abs);
            double maxV = V.Cast<double>().Max(Math.Abs);
            double gravityWaveSpeed = Math.Sqrt(g * estuaryDepth);
            double cflDt = courantNumber * Math.Min(
                dx / (maxU + gravityWaveSpeed + 1e-10),
                dy / (maxV + gravityWaveSpeed + 1e-10));
            cflDt = Math.Min(cflDt, Math.Min(
                dx * dx / (2 * (nu + eddyViscosity)),
                dy * dy / (2 * (nu + eddyViscosity))));
            return Math.Max(0.01, cflDt);
        }

        private double ApplyFluxLimiter(double q, double qUpwind, double qDownwind)
        {
            double r = (q - qUpwind + 1e-10) / (qDownwind - q + 1e-10);
            double phi = Math.Max(0, Math.Min(1, r)); // Minmod limiter
            return phi * (qDownwind - q);
        }

        public void Update(double currentTime, WetAndDryAlgo wetDryAlgo = null)
        {
            double[,] uNew = new double[nx, ny];
            double[,] vNew = new double[nx, ny];
            double[,] etaNew = new double[nx, ny];
            double[,] salinityNew = new double[nx, ny];

            double dt = ComputeCFLTimeStep();
            double tidalForcing = tidalAmplitude * Math.Sin(2 * Math.PI * currentTime / tidalPeriod);
            double tidalVelocity = tidalAmplitude * (2 * Math.PI / tidalPeriod) * Math.Cos(2 * Math.PI * currentTime / tidalPeriod);
            var (tauX, tauY) = windForcing.ComputeWindStress(); // Use WindForcing for wind stress
            double surgeLevel = stormSurgeAmplitude * Math.Exp(-currentTime / 86400.0);
            double tideSurgeInteraction = 0.05 * tidalAmplitude * stormSurgeAmplitude * Math.Sin(2 * Math.PI * currentTime / tidalPeriod);
            double waveNumber = 2 * Math.PI / (estuaryLength / 10);
            double waveVelocityX = 0.5 * waveHeight * waveNumber * Math.Sqrt(g * estuaryDepth);

            // Boundary conditions
            for (int j = 0; j < ny; j++)
            {
                uNew[0, j] = riverDischarge / (estuaryWidth * estuaryDepth);
                vNew[0, j] = 0.0;
                etaNew[0, j] = 0.0;
                salinityNew[0, j] = 0.0;
                uNew[nx - 1, j] = tidalVelocity;
                vNew[nx - 1, j] = 0.0;
                etaNew[nx - 1, j] = tidalForcing + surgeLevel + tideSurgeInteraction;
                salinityNew[nx - 1, j] = 35.0 + seasonalSalinityAmplitude * Math.Sin(2 * Math.PI * currentTime / (365 * 86400));
            }

            for (int i = 0; i < nx; i++)
            {
                vNew[i, 0] = 0.0;
                vNew[i, ny - 1] = 0.0;
                uNew[i, 0] = uNew[i, 1];
                uNew[i, ny - 1] = uNew[i, ny - 2];
                etaNew[i, 0] = etaNew[i, 1];
                etaNew[i, ny - 1] = etaNew[i, ny - 2];
                salinityNew[i, 0] = salinityNew[i, 1];
                salinityNew[i, ny - 1] = salinityNew[i, ny - 2];
            }

            // Update interior points
            for (int i = 1; i < nx - 1; i++)
                for (int j = 1; j < ny - 1; j++)
                {
                    // Skip dry cells if wetting and drying is enabled
                    if (enableWetDry && wetDryAlgo != null && !wetDryAlgo.GetWetDryStatus()[i, j])
                    {
                        uNew[i, j] = 0.0;
                        vNew[i, j] = 0.0;
                        etaNew[i, j] = Eta[i, j]; // Maintain current water level
                        salinityNew[i, j] = Salinity[i, j]; // Maintain current salinity
                        continue;
                    }

                    // Density calculation using rhoOcean
                    double rho = rho0 + Salinity[i, j] * (rhoOcean - rho0) / 35.0;

                    // Advection terms with flux limiters
                    double uAdvX = U[i, j] >= 0 ? U[i, j] * (U[i, j] - U[i - 1, j]) / dx : U[i, j] * (U[i + 1, j] - U[i, j]) / dx;
                    double uAdvY = V[i, j] >= 0 ? V[i, j] * (U[i, j] - U[i, j - 1]) / dy : V[i, j] * (U[i, j + 1] - U[i, j]) / dy;
                    double vAdvX = U[i, j] >= 0 ? U[i, j] * (V[i, j] - V[i - 1, j]) / dx : U[i, j] * (V[i + 1, j] - V[i, j]) / dx;
                    double vAdvY = V[i, j] >= 0 ? V[i, j] * (V[i, j] - V[i, j - 1]) / dy : V[i, j] * (V[i, j + 1] - V[i, j]) / dy;

                    uAdvX += U[i, j] >= 0 ? ApplyFluxLimiter(U[i, j], U[i - 1, j], U[i + 1, j]) / dx : ApplyFluxLimiter(U[i, j], U[i + 1, j], U[i - 1, j]) / dx;
                    uAdvY += V[i, j] >= 0 ? ApplyFluxLimiter(U[i, j], U[i, j - 1], U[i, j + 1]) / dy : ApplyFluxLimiter(U[i, j], U[i, j + 1], U[i, j - 1]) / dy;
                    vAdvX += U[i, j] >= 0 ? ApplyFluxLimiter(V[i, j], V[i - 1, j], V[i + 1, j]) / dx : ApplyFluxLimiter(V[i, j], V[i + 1, j], V[i - 1, j]) / dx;
                    vAdvY += V[i, j] >= 0 ? ApplyFluxLimiter(V[i, j], V[i, j - 1], V[i, j + 1]) / dy : ApplyFluxLimiter(V[i, j], V[i, j + 1], V[i, j - 1]) / dy;

                    // Diffusion terms
                    double uDiff = (nu + eddyViscosity) * (
                        (U[i + 1, j] - 2 * U[i, j] + U[i - 1, j]) / (dx * dx) +
                        (U[i, j + 1] - 2 * U[i, j] + U[i, j - 1]) / (dy * dy));
                    double vDiff = (nu + eddyViscosity) * (
                        (V[i + 1, j] - 2 * V[i, j] + V[i - 1, j]) / (dx * dx) +
                        (V[i, j + 1] - 2 * V[i, j] + V[i, j - 1]) / (dy * dy));

                    // Pressure and baroclinic terms
                    double dedx = (Eta[i + 1, j] - Eta[i - 1, j]) / (2 * dx);
                    double dedy = (Eta[i, j + 1] - Eta[i, j - 1]) / (2 * dy);
                    double drhodx = i < nx - 1 && i > 0 ?
                        (Salinity[i + 1, j] - Salinity[i - 1, j]) / (2 * dx) * (rhoOcean - rho0) / 35.0 : 0;
                    double drhody = j < ny - 1 && j > 0 ?
                        (Salinity[i, j + 1] - Salinity[i, j - 1]) / (2 * dy) * (rhoOcean - rho0) / 35.0 : 0;
                    double baroclinicX = -g * drhodx / rho0;
                    double baroclinicY = -g * drhody / rho0;

                    // Additional forcing terms
                    double coriolisU = coriolisParameter * V[i, j];
                    double coriolisV = -coriolisParameter * U[i, j];
                    double quadraticDragX = -frictionCoefficient * Math.Abs(U[i, j]) * U[i, j] / estuaryDepth;
                    double quadraticDragY = -frictionCoefficient * Math.Abs(V[i, j]) * V[i, j] / estuaryDepth;
                    double pressureGradientX = atmosphericPressureGradient / rho0;

                    // Momentum equations
                    uNew[i, j] = U[i, j] + dt * (
                        -uAdvX - uAdvY - g * dedx + uDiff + baroclinicX + coriolisU +
                        tauX / (rho0 * estuaryDepth) + waveVelocityX + quadraticDragX + pressureGradientX);
                    vNew[i, j] = V[i, j] + dt * (
                        -vAdvX - vAdvY - g * dedy + vDiff + baroclinicY + coriolisV +
                        tauY / (rho0 * estuaryDepth) + quadraticDragY);

                    // Continuity equation
                    double dudx = (U[i + 1, j] - U[i - 1, j]) / (2 * dx);
                    double dvdy = (V[i, j + 1] - V[i, j - 1]) / (2 * dy);
                    etaNew[i, j] = Eta[i, j] - dt * (estuaryDepth + Eta[i, j]) * (dudx + dvdy);

                    // Salinity transport
                    double sAdvX = U[i, j] >= 0 ? U[i, j] * (Salinity[i, j] - Salinity[i - 1, j]) / dx :
                        U[i, j] * (Salinity[i + 1, j] - Salinity[i, j]) / dx;
                    double sAdvY = V[i, j] >= 0 ? V[i, j] * (Salinity[i, j] - Salinity[i, j - 1]) / dy :
                        V[i, j] * (Salinity[i, j + 1] - Salinity[i, j]) / dy;
                    sAdvX += U[i, j] >= 0 ? ApplyFluxLimiter(Salinity[i, j], Salinity[i - 1, j], Salinity[i + 1, j]) / dx :
                        ApplyFluxLimiter(Salinity[i, j], Salinity[i + 1, j], Salinity[i - 1, j]) / dx;
                    sAdvY += V[i, j] >= 0 ? ApplyFluxLimiter(Salinity[i, j], Salinity[i, j - 1], Salinity[i, j + 1]) / dy :
                        ApplyFluxLimiter(Salinity[i, j], Salinity[i, j + 1], Salinity[i, j - 1]) / dy;
                    double sDiff = (kappa + eddyViscosity) * (
                        (Salinity[i + 1, j] - 2 * Salinity[i, j] + Salinity[i - 1, j]) / (dx * dx) +
                        (Salinity[i, j + 1] - 2 * Salinity[i, j] + Salinity[i, j - 1]) / (dy * dy));
                    double dudxShear = i < nx - 1 && i > 0 ? (U[i + 1, j] - U[i - 1, j]) / (2 * dx) : 0;
                    double dvdyShear = j < ny - 1 && j > 0 ? (V[i, j + 1] - V[i, j - 1]) / (2 * dy) : 0;
                    double shearMixing = 0.1 * Math.Abs(dudxShear + dvdyShear) * Salinity[i, j];
                    salinityNew[i, j] = Salinity[i, j] + dt * (-sAdvX - sAdvY + sDiff + shearMixing);
                    salinityNew[i, j] = Math.Max(0.0, Math.Min(35.0, salinityNew[i, j]));
                }

            // Update arrays
            U = uNew;
            V = vNew;
            Eta = etaNew;
            Salinity = salinityNew;

            // Apply wetting and drying if enabled
            if (enableWetDry && wetDryAlgo != null)
            {
                wetDryAlgo.ApplyWetDry(U, V, Eta, Salinity, dt);
            }

            // Numerical stability check
            if (U.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)) ||
                V.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)) ||
                Eta.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)) ||
                Salinity.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)))
            {
                throw new Exception("Numerical instability detected in 2D solver.");
            }
        }
    }
}