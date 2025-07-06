using System;

namespace EstuarineCirculationModeling
{
    public class EstuarineModel
    {
        public double RiverInflow { get; set; } = 0.1; // m³/s
        public double TidalAmplitude { get; set; } = 1.0; // m
        public double TidalPeriod { get; set; } = 43200; // s (12 hours)
        public double SalinityOcean { get; set; } = 35.0; // PSU
        public double TemperatureOcean { get; set; } = 20.0; // °C
        public double EstuaryLength { get; set; } = 10000; // m
        public double EstuaryDepth { get; set; } = 10.0; // m
        public double CurrentTime { get; private set; } = 0.0;
        public double SaltWedgePosition { get; private set; } = 0.0;
        public bool UseRANSSolver { get; set; } = false;
        public bool UseNonHydrostaticOverride { get; set; } = false;

        private readonly double[] salinityProfile;
        private readonly double[] temperatureProfile;
        private readonly int gridPoints = 100;
        private readonly double dt = 10.0; // Time step in seconds
        private readonly HydrodynamicSolver hydrodynamicSolver;

        public EstuarineModel()
        {
            salinityProfile = new double[gridPoints];
            temperatureProfile = new double[gridPoints];
            hydrodynamicSolver = new HydrodynamicSolver(EstuaryLength, EstuaryDepth, gridPoints, dt);
            Reset();
        }

        public void Reset()
        {
            CurrentTime = 0.0;
            SaltWedgePosition = 0.0;
            for (int i = 0; i < gridPoints; i++)
            {
                salinityProfile[i] = 0.0;
                temperatureProfile[i] = TemperatureOcean; // Initialize to ocean temperature
            }
            hydrodynamicSolver.UseNonHydrostaticOverride = UseNonHydrostaticOverride;
        }

        public void Update()
        {
            CurrentTime += dt;
            double tidalVelocity = CalculateTidalVelocity();
            double waterLevel = CalculateWaterLevel();
            double effectiveVelocity;

            if (UseRANSSolver)
            {
                hydrodynamicSolver.UseNonHydrostaticOverride = UseNonHydrostaticOverride;
                hydrodynamicSolver.Solve(RiverInflow, tidalVelocity, waterLevel, salinityProfile, temperatureProfile);
                effectiveVelocity = hydrodynamicSolver.GetVelocityAtPoint(SaltWedgePosition);
            }
            else
            {
                effectiveVelocity = tidalVelocity - (RiverInflow / (EstuaryDepth * EstuaryLength / gridPoints));
            }

            UpdateSaltWedge(effectiveVelocity);
            UpdateSalinityProfile(effectiveVelocity);
            UpdateTemperatureProfile(effectiveVelocity);
        }

        private double CalculateTidalVelocity()
        {
            return TidalAmplitude * Math.Cos(2 * Math.PI * CurrentTime / TidalPeriod);
        }

        private double CalculateWaterLevel()
        {
            return TidalAmplitude * Math.Sin(2 * Math.PI * CurrentTime / TidalPeriod);
        }

        private void UpdateSaltWedge(double velocity)
        {
            SaltWedgePosition += velocity * dt;
            SaltWedgePosition = Math.Max(0, Math.Min(SaltWedgePosition, EstuaryLength));
        }

        private void UpdateSalinityProfile(double velocity)
        {
            double dx = EstuaryLength / gridPoints;
            double diffusionCoefficient = UseRANSSolver ? hydrodynamicSolver.GetEddyViscosityAtPoint(SaltWedgePosition) : 0.1;
            double[] newSalinity = new double[gridPoints];

            for (int i = 0; i < gridPoints; i++)
            {
                double x = i * dx;
                if (x < SaltWedgePosition)
                {
                    newSalinity[i] = 0.0;
                }
                else
                {
                    newSalinity[i] = SalinityOcean;
                }

                if (i > 0 && i < gridPoints - 1)
                {
                    double advection = velocity * (salinityProfile[i + 1] - salinityProfile[i - 1]) / (2 * dx);
                    double diffusion = diffusionCoefficient * (salinityProfile[i + 1] - 2 * salinityProfile[i] + salinityProfile[i - 1]) / (dx * dx);
                    newSalinity[i] -= dt * advection;
                    newSalinity[i] += dt * diffusion;
                }
            }

            for (int i = 0; i < gridPoints; i++)
            {
                salinityProfile[i] = Math.Max(0, Math.Min(SalinityOcean, newSalinity[i]));
            }
        }

        private void UpdateTemperatureProfile(double velocity)
        {
            double dx = EstuaryLength / gridPoints;
            double diffusionCoefficient = UseRANSSolver ? hydrodynamicSolver.GetEddyViscosityAtPoint(SaltWedgePosition) : 0.1;
            double[] newTemperature = new double[gridPoints];

            for (int i = 0; i < gridPoints; i++)
            {
                double x = i * dx;
                if (x < SaltWedgePosition)
                {
                    newTemperature[i] = 15.0; // River temperature
                }
                else
                {
                    newTemperature[i] = TemperatureOcean;
                }

                if (i > 0 && i < gridPoints - 1)
                {
                    double advection = velocity * (temperatureProfile[i + 1] - temperatureProfile[i - 1]) / (2 * dx);
                    double diffusion = diffusionCoefficient * (temperatureProfile[i + 1] - 2 * temperatureProfile[i] + temperatureProfile[i - 1]) / (dx * dx);
                    newTemperature[i] = temperatureProfile[i] - dt * advection + dt * diffusion;
                }
            }

            for (int i = 0; i < gridPoints; i++)
            {
                temperatureProfile[i] = Math.Max(0, Math.Min(TemperatureOcean, newTemperature[i]));
            }
        }

        public double GetSalinityAtPoint(double x)
        {
            int index = (int)(x / EstuaryLength * gridPoints);
            index = Math.Max(0, Math.Min(gridPoints - 1, index));
            return salinityProfile[index];
        }

        public double GetTemperatureAtPoint(double x)
        {
            int index = (int)(x / EstuaryLength * gridPoints);
            index = Math.Max(0, Math.Min(gridPoints - 1, index));
            return temperatureProfile[index];
        }

        public double GetMaxSalinity()
        {
            double maxSalinity = 0.0;
            for (int i = 0; i < gridPoints; i++)
            {
                if (salinityProfile[i] > maxSalinity)
                {
                    maxSalinity = salinityProfile[i];
                }
            }
            return maxSalinity;
        }

        public double[] GetVelocityProfile()
        {
            return hydrodynamicSolver.GetVelocityProfile();
        }

        public double[] GetSalinityProfile() // New method to access salinity profile
        {
            return (double[])salinityProfile.Clone();
        }

        public double[] GetTemperatureProfile() // New method to access temperature profile
        {
            return (double[])temperatureProfile.Clone();
        }

        public double GetEddyViscosityAtPoint(double x)
        {
            return hydrodynamicSolver.GetEddyViscosityAtPoint(x);
        }

        public bool IsNonHydrostatic()
        {
            return hydrodynamicSolver.IsNonHydrostatic;
        }
    }
}