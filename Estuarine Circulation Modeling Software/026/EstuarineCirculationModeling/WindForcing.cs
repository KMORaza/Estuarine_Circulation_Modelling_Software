using System;

namespace EstuarineCirculationModeling
{
    public class WindForcing
    {
        private readonly double rhoAir = 1.225; // Air density (kg/m^3)
        private double windSpeed; // Wind speed at 10m (m/s)
        private double windDirection; // Wind direction (degrees)
        private double waveHeight; // Significant wave height (m)

        public WindForcing(double windSpeed, double windDirection, double waveHeight)
        {
            this.windSpeed = Math.Max(0.0, windSpeed);
            this.windDirection = windDirection % 360.0;
            this.waveHeight = Math.Max(0.0, waveHeight);
        }

        // Compute wind drag coefficient based on wind speed and sea state (wave height)
        public double ComputeDragCoefficient()
        {
            // Simplified parameterization: C_d = (0.75 + 0.067 * U10 + 0.1 * H_s) * 10^-3
            double Cd = (0.75 + 0.067 * windSpeed + 0.1 * waveHeight) * 1e-3;
            // Bound C_d between 0.001 and 0.003 for physical realism
            return Math.Max(0.001, Math.Min(0.003, Cd));
        }

        // Compute wind stress components (tau_x, tau_y) in N/m^2
        public (double tauX, double tauY) ComputeWindStress()
        {
            double Cd = ComputeDragCoefficient();
            double windSpeedSquared = windSpeed * windSpeed;
            double tauX = rhoAir * Cd * windSpeedSquared * Math.Cos(windDirection * Math.PI / 180.0);
            double tauY = rhoAir * Cd * windSpeedSquared * Math.Sin(windDirection * Math.PI / 180.0);
            return (tauX, tauY);
        }

        // Update wind parameters (called when UI inputs change)
        public void UpdateParameters(double windSpeed, double windDirection, double waveHeight)
        {
            this.windSpeed = Math.Max(0.0, windSpeed);
            this.windDirection = windDirection % 360.0;
            this.waveHeight = Math.Max(0.0, waveHeight);
        }
    }
}