using System;
using System.Collections.Generic;

namespace EstuarineCirculationModeling
{
    public class SimpsonHunterMechanismParam
    {
        private readonly double gravity = 9.81; // m/s²
        private readonly double referenceDensity = 1000.0; // kg/m³
        private readonly int numSigmaLayers;
        private readonly double waveAmplitude; // Wave amplitude for Stokes drift (m)
        private readonly double wavePeriod; // Wave period for Stokes drift (s)
        private readonly double internalTideAmplitude; // Amplitude of internal tide perturbation (m/s)
        private readonly double simpsonHunterCoefficient; // Coefficient for tidal straining effect

        public SimpsonHunterMechanismParam(int numSigmaLayers, double waveAmplitude = 0.5, double wavePeriod = 10.0, double internalTideAmplitude = 0.05, double simpsonHunterCoefficient = 0.1)
        {
            this.numSigmaLayers = numSigmaLayers;
            this.waveAmplitude = waveAmplitude; // Default wave amplitude
            this.wavePeriod = wavePeriod; // Default wave period
            this.internalTideAmplitude = internalTideAmplitude; // Default internal tide amplitude
            this.simpsonHunterCoefficient = simpsonHunterCoefficient; // Default tidal straining coefficient
        }

        // Compute Stokes drift velocity (x-direction) at a given depth and tidal phase
        public double ComputeStokesDrift(double depth, double sigma, double tidalPhase)
        {
            // Stokes drift: u_s = (a * omega)^2 / (k * H) * cosh(2k(H-z)) / sinh^2(kH)
            // Simplified for shallow water: u_s = a * omega * exp(-2kz)
            double k = 2 * Math.PI / (wavePeriod * Math.Sqrt(gravity * depth)); // Wave number
            double z = sigma * depth; // Depth coordinate
            double stokesDrift = waveAmplitude * 2 * Math.PI / wavePeriod * Math.Exp(-2 * k * z) * Math.Sin(tidalPhase);
            return Math.Max(-0.5, Math.Min(0.5, stokesDrift)); // Cap for stability
        }

        // Compute internal tide effect on vertical velocity and TKE production
        public double ComputeInternalTideEffect(AsymmTidalMix.Cell cell, int k, double depth, double sigmaStep, double tidalPhase, out double tkeProduction)
        {
            tkeProduction = 0.0;
            if (k >= numSigmaLayers - 1) return 0.0; // No effect at top layer

            // Internal tide: vertical velocity perturbation due to stratification and bathymetry
            double sigma = k * sigmaStep + sigmaStep / 2;
            double z = sigma * depth;
            double dRho_dz = 0.0;
            if (k < numSigmaLayers - 1)
            {
                double rho1 = referenceDensity + 0.8 * cell.Salinity[k];
                double rho2 = referenceDensity + 0.8 * cell.Salinity[k + 1];
                dRho_dz = (rho2 - rho1) / (sigmaStep * depth);
            }
            double N2 = -gravity / referenceDensity * dRho_dz;
            N2 = Math.Max(0.0, Math.Min(1e-3, N2)); // Cap buoyancy frequency squared

            // Internal tide vertical velocity perturbation
            double wTide = internalTideAmplitude * Math.Sin(tidalPhase) * Math.Sqrt(N2) * Math.Cos(2 * Math.PI * z / depth);
            wTide = Math.Max(-0.2, Math.Min(0.2, wTide)); // Cap for stability

            // Additional TKE production due to internal tide shear
            if (k < numSigmaLayers - 1)
            {
                double dw_dz = (cell.W[k] - cell.W[k + 1]) / (sigmaStep * depth);
                tkeProduction = 0.1 * internalTideAmplitude * dw_dz * dw_dz; // Simplified shear production
                tkeProduction = Math.Max(0.0, Math.Min(1e-4, tkeProduction)); // Cap for stability
            }

            return wTide;
        }

        // Compute tidal straining (Simpson-Hunter) effect on salinity gradient
        public double ComputeTidalStraining(AsymmTidalMix.Cell cell, int k, double depth, double sigmaStep, double tidalVelocity, List<AsymmTidalMix.Cell> cells, int cellIdx, double avgCellWidth)
        {
            if (k >= numSigmaLayers - 1 || cellIdx >= cells.Count - 1) return 0.0; // No effect at top layer or boundary

            // Tidal straining: u * ds/dx enhances or suppresses stratification
            double ds_dx = 0.0;
            if (cellIdx > 0 && cellIdx < cells.Count - 1)
            {
                ds_dx = (cells[cellIdx + 1].Salinity[k] - cells[cellIdx - 1].Salinity[k]) / (2 * avgCellWidth);
            }
            else if (cellIdx == 0)
            {
                ds_dx = (cells[1].Salinity[k] - cell.Salinity[k]) / avgCellWidth;
            }
            else
            {
                ds_dx = (cell.Salinity[k] - cells[cellIdx - 1].Salinity[k]) / avgCellWidth;
            }
            ds_dx = Math.Max(-10.0, Math.Min(10.0, ds_dx)); // Cap salinity gradient

            // Straining term: -u * ds/dx
            double straining = -simpsonHunterCoefficient * tidalVelocity * ds_dx;
            straining = Math.Max(-1.0, Math.Min(1.0, straining)); // Cap for stability

            return straining;
        }
    }
}