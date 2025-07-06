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
    }
}