using System;

namespace EstuarineCirculationModeling
{
    public class WaveCurrentInteraction
    {
        private double waveHeight; // Wave height (m)
        private double waveDirection; // Wave direction (degrees, aligned with wind for simplicity)
        private double wavePeriod; // Wave period (s)
        private double depth; // Water depth (m)
        private readonly double g = 9.81; // Gravitational acceleration (m/s²)

        public WaveCurrentInteraction(double waveHeight, double waveDirection, double wavePeriod, double depth)
        {
            this.waveHeight = waveHeight;
            this.waveDirection = waveDirection;
            this.wavePeriod = wavePeriod;
            this.depth = depth;
        }

        // Update parameters if wave conditions change
        public void UpdateParameters(double waveHeight, double waveDirection, double wavePeriod, double depth)
        {
            this.waveHeight = waveHeight;
            this.waveDirection = waveDirection;
            this.wavePeriod = wavePeriod;
            this.depth = depth;
        }

        // Compute Stokes drift velocities (uStokes, vStokes) based on linear wave theory
        public (double uStokes, double vStokes) ComputeStokesDrift()
        {
            if (waveHeight <= 0 || wavePeriod <= 0) return (0.0, 0.0);

            // Wave number k = 2π / wavelength, using shallow water dispersion relation
            double waveLength = Math.Sqrt(g * depth) * wavePeriod; // Approximate for shallow water
            double k = 2 * Math.PI / waveLength;
            double omega = 2 * Math.PI / wavePeriod; // Angular frequency
            double a = waveHeight / 2; // Wave amplitude

            // Stokes drift velocity: u_s = a²ωk * cosh(2kz) / (2 * sinh²(kh)) ≈ a²ωk for shallow water (z ≈ 0)
            double uStokesMag = (a * a * omega * k) / (2 * Math.Sinh(k * depth) * Math.Sinh(k * depth)) * Math.Cosh(2 * k * 0); // Near surface
            if (double.IsNaN(uStokesMag) || double.IsInfinity(uStokesMag)) uStokesMag = 0.0;

            // Resolve into x and y components based on wave direction
            double dirRad = waveDirection * Math.PI / 180.0;
            double uStokes = uStokesMag * Math.Cos(dirRad);
            double vStokes = uStokesMag * Math.Sin(dirRad);

            return (uStokes, vStokes);
        }

        // Compute wave-enhanced bottom friction coefficient (Grant-Madsen simplified)
        public double ComputeWaveEnhancedFriction(double baseFrictionCoefficient)
        {
            if (waveHeight <= 0 || wavePeriod <= 0) return baseFrictionCoefficient;

            // Wave orbital velocity at bottom: U_b = aω / sinh(kh)
            double waveLength = Math.Sqrt(g * depth) * wavePeriod;
            double k = 2 * Math.PI / waveLength;
            double omega = 2 * Math.PI / wavePeriod;
            double a = waveHeight / 2;
            double ub = a * omega / Math.Sinh(k * depth);
            if (double.IsNaN(ub) || double.IsInfinity(ub)) ub = 0.0;

            // Wave-enhanced friction coefficient: C_d = C_d0 * (1 + β * U_b), β = empirical constant
            const double beta = 0.2; // Empirical factor for wave-current interaction
            double enhancedFriction = baseFrictionCoefficient * (1.0 + beta * Math.Abs(ub));
            return Math.Max(baseFrictionCoefficient, Math.Min(0.01, enhancedFriction)); // Cap to avoid instability
        }
    }
}