using System;

namespace EstuarineCirculationModeling
{
    public class HydrodynamicSolver
    {
        private readonly int gridPoints;
        private readonly double dx; // Spatial step
        private readonly double dt; // Time step
        private readonly double estuaryLength;
        private readonly double estuaryDepth;
        private readonly double[] velocity; // Velocity field (m/s)
        private readonly double[] pressure; // Pressure field (Pa)
        private readonly double[] k; // Turbulent kinetic energy (m²/s²)
        private readonly double[] epsilon; // Dissipation rate (m²/s³)
        private readonly double[] eddyViscosity; // Turbulent eddy viscosity (m²/s)
        private readonly double kinematicViscosity = 1e-6; // Kinematic viscosity (m²/s)
        private readonly double coriolisParameter = 1e-4; // Coriolis parameter (s^-1)
        private readonly double density = 1000.0; // Reference density (kg/m³)
        private readonly double g = 9.81; // Gravitational acceleration (m/s²)
        private readonly double resolutionThreshold = 200.0; // dx threshold (m) for non-hydrostatic
        private readonly BaroclinicFlow baroclinicFlow; // Baroclinic flow calculator
        private bool useNonHydrostatic; // Manual override for approximation

        // k-ε model constants
        private readonly double Cmu = 0.09;
        private readonly double sigmaK = 1.0;
        private readonly double sigmaEpsilon = 1.3;
        private readonly double C1e = 1.44;
        private readonly double C2e = 1.92;

        public bool IsNonHydrostatic { get; private set; }
        public bool UseNonHydrostaticOverride { get; set; }

        public HydrodynamicSolver(double estuaryLength, double estuaryDepth, int gridPoints, double dt)
        {
            this.estuaryLength = estuaryLength;
            this.estuaryDepth = estuaryDepth;
            this.gridPoints = gridPoints;
            this.dx = estuaryLength / gridPoints;
            this.dt = dt;
            velocity = new double[gridPoints];
            pressure = new double[gridPoints];
            k = new double[gridPoints];
            epsilon = new double[gridPoints];
            eddyViscosity = new double[gridPoints];
            useNonHydrostatic = dx <= resolutionThreshold; // Default: non-hydrostatic if dx ≤ 200m
            baroclinicFlow = new BaroclinicFlow(estuaryLength, estuaryDepth, gridPoints);
            InitializeFields();
        }

        private void InitializeFields()
        {
            for (int i = 0; i < gridPoints; i++)
            {
                velocity[i] = 0.0;
                pressure[i] = density * g * estuaryDepth; // Initial hydrostatic pressure
                k[i] = 1e-3;
                epsilon[i] = 1e-6;
                eddyViscosity[i] = 0.01;
            }
        }

        public void Solve(double riverInflow, double tidalVelocity, double waterLevel, double[] salinity, double[] temperature)
        {
            double[] newVelocity = new double[gridPoints];
            double[] newK = new double[gridPoints];
            double[] newEpsilon = new double[gridPoints];
            IsNonHydrostatic = UseNonHydrostaticOverride || dx <= resolutionThreshold;

            // Compute baroclinic pressure gradient
            double[] baroclinicGradient = baroclinicFlow.ComputeBaroclinicGradient(salinity, temperature);

            // Step 1: Solve momentum equation
            for (int i = 1; i < gridPoints - 1; i++)
            {
                double advection = velocity[i] * (velocity[i + 1] - velocity[i - 1]) / (2 * dx);
                double pressureGradient = -(pressure[i + 1] - pressure[i - 1]) / (2 * dx * density);
                double baroclinicTerm = -baroclinicGradient[i] / density; // Add baroclinic term
                double totalViscosity = kinematicViscosity + eddyViscosity[i];
                double diffusion = totalViscosity * (velocity[i + 1] - 2 * velocity[i] + velocity[i - 1]) / (dx * dx);
                double coriolis = coriolisParameter * velocity[i];

                newVelocity[i] = velocity[i] + dt * (-advection + pressureGradient + baroclinicTerm + diffusion - coriolis);
            }

            // Boundary conditions
            newVelocity[0] = riverInflow / estuaryDepth;
            newVelocity[gridPoints - 1] = tidalVelocity;

            // Step 2: Solve k-ε turbulence model
            for (int i = 1; i < gridPoints - 1; i++)
            {
                double shear = Math.Abs((velocity[i + 1] - velocity[i - 1]) / (2 * dx));
                double production = eddyViscosity[i] * shear * shear;

                double kAdvection = velocity[i] * (k[i + 1] - k[i - 1]) / (2 * dx);
                double kDiffusion = (kinematicViscosity + eddyViscosity[i] / sigmaK) * (k[i + 1] - 2 * k[i] + k[i - 1]) / (dx * dx);
                newK[i] = k[i] + dt * (-kAdvection + production - epsilon[i] + kDiffusion);

                double epsilonAdvection = velocity[i] * (epsilon[i + 1] - epsilon[i - 1]) / (2 * dx);
                double epsilonDiffusion = (kinematicViscosity + eddyViscosity[i] / sigmaEpsilon) * (epsilon[i + 1] - 2 * epsilon[i] + epsilon[i - 1]) / (dx * dx);
                newEpsilon[i] = epsilon[i] + dt * (-epsilonAdvection + C1e * (epsilon[i] / k[i]) * production - C2e * (epsilon[i] * epsilon[i] / k[i]) + epsilonDiffusion);

                eddyViscosity[i] = Cmu * newK[i] * newK[i] / Math.Max(newEpsilon[i], 1e-10);
            }

            newK[0] = newK[1];
            newK[gridPoints - 1] = newK[gridPoints - 2];
            newEpsilon[0] = newEpsilon[1];
            newEpsilon[gridPoints - 1] = newEpsilon[gridPoints - 2];

            // Step 3: Pressure update
            double[] newPressure = new double[gridPoints];
            if (IsNonHydrostatic)
            {
                for (int i = 1; i < gridPoints - 1; i++)
                {
                    double divergence = (newVelocity[i + 1] - newVelocity[i - 1]) / (2 * dx);
                    newPressure[i] = pressure[i] - density * divergence * dx * dx / dt;
                }
                newPressure[0] = newPressure[1];
                newPressure[gridPoints - 1] = density * g * (estuaryDepth + waterLevel);
            }
            else
            {
                for (int i = 0; i < gridPoints; i++)
                {
                    newPressure[i] = density * g * (estuaryDepth + waterLevel * (i * dx / estuaryLength));
                }
            }

            Array.Copy(newVelocity, velocity, gridPoints);
            Array.Copy(newPressure, pressure, gridPoints);
            Array.Copy(newK, k, gridPoints);
            Array.Copy(newEpsilon, epsilon, gridPoints);
        }

        public double GetVelocityAtPoint(double x)
        {
            int index = (int)(x / estuaryLength * gridPoints);
            index = Math.Max(0, Math.Min(gridPoints - 1, index));
            return velocity[index];
        }

        public double[] GetVelocityProfile()
        {
            return (double[])velocity.Clone();
        }

        public double GetEddyViscosityAtPoint(double x)
        {
            int index = (int)(x / estuaryLength * gridPoints);
            index = Math.Max(0, Math.Min(gridPoints - 1, index));
            return eddyViscosity[index];
        }
    }
}