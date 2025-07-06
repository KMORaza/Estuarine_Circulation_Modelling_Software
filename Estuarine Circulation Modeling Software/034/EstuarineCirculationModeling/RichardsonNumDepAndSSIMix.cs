using System;
using System.Collections.Generic;

namespace EstuarineCirculationModeling
{
    public class RichardsonNumDepAndSSIMix
    {
        private readonly double gravity = 9.81; // m/s²
        private readonly double referenceDensity = 1000.0; // kg/m³
        private readonly double kinematicViscosity = 1e-6; // m²/s
        private readonly int numSigmaLayers;
        private readonly double cMu = 0.09; // k-ε model constant

        public RichardsonNumDepAndSSIMix(int numSigmaLayers)
        {
            this.numSigmaLayers = numSigmaLayers;
        }

        // Compute Richardson number and adjust turbulent viscosity
        public double ComputeAdjustedTurbulentViscosity(AsymmTidalMix.Cell cell, int k, double depth, double sigmaStep, out double richardsonNumber)
        {
            richardsonNumber = 0.0;
            double nuT = cMu * cell.K[k] * cell.K[k] / cell.Epsilon[k];

            // Avoid division by zero or invalid values
            if (double.IsNaN(nuT) || double.IsInfinity(nuT) || cell.Epsilon[k] < 1e-8 || cell.K[k] < 1e-6)
            {
                nuT = kinematicViscosity;
                return Math.Max(1e-6, Math.Min(1e-2, nuT));
            }

            // Compute buoyancy frequency squared (N²)
            double dRho_dz = 0.0;
            if (k < numSigmaLayers - 1)
            {
                double rho1 = referenceDensity + 0.8 * cell.Salinity[k];
                double rho2 = referenceDensity + 0.8 * cell.Salinity[k + 1];
                dRho_dz = (rho2 - rho1) / (sigmaStep * depth);
            }
            double N2 = -gravity / referenceDensity * dRho_dz;
            N2 = Math.Max(0.0, Math.Min(1e-3, N2)); // Cap N² to prevent instability

            // Compute shear squared (S²)
            double S2 = cell.Shear[k] * cell.Shear[k];
            S2 = Math.Max(1e-6, Math.Min(100.0, S2)); // Cap shear squared

            // Compute Richardson number
            richardsonNumber = N2 / S2;
            richardsonNumber = double.IsNaN(richardsonNumber) || double.IsInfinity(richardsonNumber) ? 0.0 : Math.Max(0.0, Math.Min(10.0, richardsonNumber));

            // Stability function: reduces mixing when Ri > 0.25
            double fRi = 1.0 / (1.0 + 10.0 * richardsonNumber);
            fRi = Math.Max(0.1, Math.Min(1.0, fRi)); // Ensure f(Ri) is reasonable
            nuT *= fRi;

            // Cap turbulent viscosity
            return Math.Max(1e-6, Math.Min(1e-2, nuT));
        }

        // Compute shear-strain-induced TKE production
        public double ComputeShearStrainProduction(AsymmTidalMix.Cell cell, int k, double depth, double sigmaStep, double nuT)
        {
            if (k >= numSigmaLayers - 1) return 0.0; // No shear at top layer

            // Compute strain tensor components (simplified to vertical shear)
            double du_dz = (cell.U[k] - cell.U[k + 1]) / (sigmaStep * depth);
            double dv_dz = (cell.V[k] - cell.V[k + 1]) / (sigmaStep * depth);
            du_dz = Math.Max(-100.0, Math.Min(100.0, du_dz));
            dv_dz = Math.Max(-100.0, Math.Min(100.0, dv_dz));

            // Shear-strain production: P = nuT * (du_i/dx_j + du_j/dx_i) * du_i/dx_j
            // Simplified to vertical shear terms only
            double shearStrainProd = nuT * (du_dz * du_dz + dv_dz * dv_dz);
            shearStrainProd = Math.Max(0.0, Math.Min(1e-3, shearStrainProd)); // Cap to prevent instability

            return shearStrainProd;
        }

        // Compute average Richardson number across cells for output
        public double ComputeAverageRichardsonNumber(List<AsymmTidalMix.Cell> cells, double avgCellWidth)
        {
            double totalRi = 0.0;
            int validCells = 0;

            foreach (var cell in cells)
            {
                double depth = cell.Depth;
                double sigmaStep = 1.0 / numSigmaLayers;
                for (int k = 0; k < numSigmaLayers - 1; k++)
                {
                    double dRho_dz = 0.0;
                    double rho1 = referenceDensity + 0.8 * cell.Salinity[k];
                    double rho2 = referenceDensity + 0.8 * cell.Salinity[k + 1];
                    dRho_dz = (rho2 - rho1) / (sigmaStep * depth);
                    double N2 = -gravity / referenceDensity * dRho_dz;
                    N2 = Math.Max(0.0, Math.Min(1e-3, N2));

                    double S2 = cell.Shear[k] * cell.Shear[k];
                    S2 = Math.Max(1e-6, Math.Min(100.0, S2));

                    double ri = N2 / S2;
                    if (!double.IsNaN(ri) && !double.IsInfinity(ri))
                    {
                        totalRi += Math.Max(0.0, Math.Min(10.0, ri));
                        validCells++;
                    }
                }
            }

            return validCells > 0 ? totalRi / validCells : 0.0;
        }
    }
}