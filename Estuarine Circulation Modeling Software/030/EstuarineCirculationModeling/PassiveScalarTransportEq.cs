using System;

namespace EstuarineCirculationModeling
{
    public class PassiveScalarTransportEq
    {
        private readonly int gridPoints;
        private readonly double dx; // Spatial step
        private readonly double estuaryLength;

        public PassiveScalarTransportEq(double estuaryLength, int gridPoints)
        {
            this.gridPoints = gridPoints;
            this.dx = estuaryLength / gridPoints;
            this.estuaryLength = estuaryLength;
        }

        public double[] SolvePassiveScalarTransport(double[] concentration, double[] velocity, double dt, double riverConcentration, double oceanConcentration, double saltWedgePosition, double estuaryLength, double[] eddyDiffusivity)
        {
            double[] newConcentration = new double[gridPoints];

            // Apply advection-diffusion equation using finite difference
            for (int i = 1; i < gridPoints - 1; i++)
            {
                double advection = velocity[i] * (concentration[i + 1] - concentration[i - 1]) / (2 * dx);
                double diffusion = (eddyDiffusivity[i + 1] * (concentration[i + 1] - concentration[i]) / dx - eddyDiffusivity[i] * (concentration[i] - concentration[i - 1]) / dx) / dx;
                newConcentration[i] = concentration[i] + dt * (-advection + diffusion);
            }

            // Boundary conditions
            newConcentration[0] = riverConcentration; // River boundary
            newConcentration[gridPoints - 1] = oceanConcentration; // Ocean boundary

            // Apply mixing near salt wedge
            double saltWedgeIndex = saltWedgePosition / dx;
            for (int i = 0; i < gridPoints; i++)
            {
                double distance = Math.Abs(i * dx - saltWedgePosition);
                double mixingFactor = Math.Exp(-distance * distance / (0.1 * estuaryLength * estuaryLength));
                newConcentration[i] = (1 - mixingFactor) * newConcentration[i] + mixingFactor * concentration[i];
            }

            return newConcentration;
        }
    }
}