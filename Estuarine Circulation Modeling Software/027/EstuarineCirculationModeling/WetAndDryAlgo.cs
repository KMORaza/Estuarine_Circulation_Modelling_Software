using System;

namespace EstuarineCirculationModeling
{
    public class WetAndDryAlgo
    {
        private double[,] bathymetry; // Bathymetry (positive downward)
        private bool[,] isWet; // Wet/dry status
        private double Dmin; // Minimum depth threshold
        private int nx, ny;
        private double dx, dy;
        private double estuaryLength, estuaryWidth;

        public WetAndDryAlgo(double[,] bathymetry, double Dmin, int nx, int ny, double estuaryLength, double estuaryWidth)
        {
            this.bathymetry = bathymetry;
            this.Dmin = Dmin;
            this.nx = nx;
            this.ny = ny;
            this.estuaryLength = estuaryLength;
            this.estuaryWidth = estuaryWidth;
            this.dx = estuaryLength / (nx - 1);
            this.dy = estuaryWidth / (ny - 1);
            isWet = new bool[nx, ny];
            InitializeWetDry(new double[nx, ny]); // Initialize with zero water level
        }

        public void InitializeWetDry(double[,] eta)
        {
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    double totalDepth = eta[i, j] + bathymetry[i, j]; // eta is water level, bathymetry is positive downward
                    isWet[i, j] = totalDepth >= Dmin;
                }
        }

        public void ApplyWetDry(double[,] u, double[,] v, double[,] eta, double[,] salinity, double dt)
        {
            // Update wet/dry status
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    double totalDepth = eta[i, j] + bathymetry[i, j];
                    isWet[i, j] = totalDepth >= Dmin;
                    if (!isWet[i, j])
                    {
                        u[i, j] = 0.0;
                        v[i, j] = 0.0;
                        eta[i, j] = -bathymetry[i, j]; // Set water level to bed level
                        salinity[i, j] = 0.0; // No salinity in dry cells
                    }
                }

            // Adjust fluxes at wet/dry interfaces
            for (int i = 1; i < nx - 1; i++)
                for (int j = 1; j < ny - 1; j++)
                {
                    if (!isWet[i, j]) continue;

                    // Check neighboring cells
                    if (!isWet[i + 1, j])
                    {
                        u[i, j] = Math.Min(u[i, j], 0.0); // Prevent flow into dry cell
                    }
                    if (!isWet[i - 1, j])
                    {
                        u[i - 1, j] = Math.Max(u[i - 1, j], 0.0); // Prevent flow from dry cell
                    }
                    if (!isWet[i, j + 1])
                    {
                        v[i, j] = Math.Min(v[i, j], 0.0); // Prevent flow into dry cell
                    }
                    if (!isWet[i, j - 1])
                    {
                        v[i, j - 1] = Math.Max(v[i, j - 1], 0.0); // Prevent flow from dry cell
                    }
                }

            // Mass conservation correction
            for (int i = 1; i < nx - 1; i++)
                for (int j = 1; j < ny - 1; j++)
                {
                    if (!isWet[i, j]) continue;

                    double fluxIn = 0.0, fluxOut = 0.0;
                    if (isWet[i - 1, j])
                        fluxIn += u[i - 1, j] * (eta[i - 1, j] + bathymetry[i - 1, j]) * dy * dt;
                    if (isWet[i + 1, j])
                        fluxOut += u[i, j] * (eta[i, j] + bathymetry[i, j]) * dy * dt;
                    if (isWet[i, j - 1])
                        fluxIn += v[i, j - 1] * (eta[i, j - 1] + bathymetry[i, j - 1]) * dx * dt;
                    if (isWet[i, j + 1])
                        fluxOut += v[i, j] * (eta[i, j] + bathymetry[i, j]) * dx * dt;

                    eta[i, j] += (fluxIn - fluxOut) / (dx * dy);
                    eta[i, j] = Math.Max(-bathymetry[i, j], eta[i, j]); // Ensure no negative depth
                }
        }

        public bool[,] GetWetDryStatus()
        {
            return isWet;
        }
    }
}