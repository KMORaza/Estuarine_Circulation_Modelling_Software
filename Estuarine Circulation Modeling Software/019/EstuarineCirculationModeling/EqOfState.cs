using System;

namespace EstuarineCirculationModeling
{
    public static class EqOfState
    {
        // Compute seawater density (kg/m³) using UNESCO 1983 equation of state
        public static double ComputeDensity(double salinity, double temperature, double pressure)
        {
            // Salinity (PSU), temperature (°C), pressure (dbar, 1 dbar = 10^4 Pa)
            // Ensure physical bounds
            salinity = Math.Max(0, Math.Min(40, salinity)); // PSU
            temperature = Math.Max(-2, Math.Min(40, temperature)); // °C
            pressure = Math.Max(0, pressure); // dbar

            // UNESCO 1983 coefficients for pure water density (ρ₀ at P=0)
            double a0 = 999.842594;
            double a1 = 6.793952e-2;
            double a2 = -9.095290e-3;
            double a3 = 1.001685e-4;
            double a4 = -1.120083e-6;
            double a5 = 6.536332e-9;

            // Salinity coefficients
            double b0 = 8.24493e-1;
            double b1 = -4.0899e-3;
            double b2 = 7.6438e-5;
            double b3 = -8.2467e-7;
            double b4 = 5.3875e-9;

            // Temperature-salinity interaction
            double c0 = -5.72466e-3;
            double c1 = 1.0227e-4;
            double c2 = -1.6546e-6;

            // Salinity^1.5 term
            double d0 = 4.8314e-4;

            // Pressure coefficients (bulk modulus)
            double e0 = 19652.21;
            double e1 = 148.4206;
            double e2 = -2.327105;
            double e3 = 1.360477e-2;
            double e4 = -5.155288e-5;

            double f0 = 54.6746;
            double f1 = -0.603459;
            double f2 = 1.09987e-2;
            double f3 = -6.1670e-5;

            double g0 = 7.944e-2;
            double g1 = 1.6483e-2;
            double g2 = -5.3009e-4;

            // Pure water density at P=0
            double rho0 = a0 + a1 * temperature + a2 * Math.Pow(temperature, 2) +
                          a3 * Math.Pow(temperature, 3) + a4 * Math.Pow(temperature, 4) +
                          a5 * Math.Pow(temperature, 5);

            // Salinity contribution
            double salinityTerm = b0 + b1 * temperature + b2 * Math.Pow(temperature, 2) +
                                  b3 * Math.Pow(temperature, 3) + b4 * Math.Pow(temperature, 4);

            // Temperature-salinity interaction
            double interactionTerm = c0 + c1 * temperature + c2 * Math.Pow(temperature, 2);

            // Density at atmospheric pressure
            double rho = rho0 + salinity * salinityTerm + salinity * Math.Sqrt(salinity) * d0 +
                         salinity * interactionTerm;

            // Compressibility (secant bulk modulus)
            double K0 = e0 + e1 * temperature + e2 * Math.Pow(temperature, 2) +
                        e3 * Math.Pow(temperature, 3) + e4 * Math.Pow(temperature, 4);

            double K1 = f0 + f1 * temperature + f2 * Math.Pow(temperature, 2) +
                        f3 * Math.Pow(temperature, 3);

            double K2 = g0 + g1 * temperature + g2 * Math.Pow(temperature, 2);

            double K = K0 + salinity * K1 + salinity * Math.Sqrt(salinity) * K2;

            // Density with pressure correction
            rho = rho * (1 + pressure / (K - pressure));

            return rho;
        }
    }
}