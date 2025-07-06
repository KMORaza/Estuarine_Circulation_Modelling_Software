using System;

namespace EstuarineCirculationModeling
{
    public abstract class VerticalDiscretization
    {
        protected int gridPointsX;
        protected int gridPointsZ;
        protected double dxi;
        protected double deta;
        protected double estuaryLength;
        protected double estuaryDepth;
        protected double[,] x;
        protected double[,] z;
        protected double[,] z_xi;
        protected double[,] z_eta;

        public double[,] Z => z;
        public double[,] Z_xi => z_xi;
        public double[,] Z_eta => z_eta;

        public VerticalDiscretization(int gridPointsX, int gridPointsZ, double dxi, double deta, double estuaryLength, double estuaryDepth, double[,] x)
        {
            this.gridPointsX = gridPointsX;
            this.gridPointsZ = gridPointsZ;
            this.dxi = dxi;
            this.deta = deta;
            this.estuaryLength = estuaryLength;
            this.estuaryDepth = estuaryDepth;
            this.x = x;
            z = new double[gridPointsX, gridPointsZ];
            z_xi = new double[gridPointsX, gridPointsZ];
            z_eta = new double[gridPointsX, gridPointsZ];
        }

        public abstract void InitializeGrid();
        public abstract void UpdateDepth(double newEstuaryDepth);
    }

    public class SigmaCoordinates : VerticalDiscretization
    {
        public SigmaCoordinates(int gridPointsX, int gridPointsZ, double dxi, double deta, double estuaryLength, double estuaryDepth, double[,] x)
            : base(gridPointsX, gridPointsZ, dxi, deta, estuaryLength, estuaryDepth, x)
        {
            InitializeGrid();
        }

        public override void InitializeGrid()
        {
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    double eta = j * deta;
                    z[i, j] = eta * estuaryDepth; // Sigma: z follows depth
                }
            }

            ComputeMetricTerms();
        }

        public override void UpdateDepth(double newEstuaryDepth)
        {
            estuaryDepth = newEstuaryDepth;
            InitializeGrid();
        }

        private void ComputeMetricTerms()
        {
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    z_xi[i, j] = i > 0 && i < gridPointsX - 1 ? (z[i + 1, j] - z[i - 1, j]) / (2 * dxi) :
                        i == 0 ? (z[1, j] - z[0, j]) / dxi : (z[gridPointsX - 1, j] - z[gridPointsX - 2, j]) / dxi;
                    z_eta[i, j] = j > 0 && j < gridPointsZ - 1 ? (z[i, j + 1] - z[i, j - 1]) / (2 * deta) :
                        j == 0 ? (z[i, 1] - z[i, 0]) / deta : (z[i, gridPointsZ - 1] - z[i, gridPointsZ - 2]) / deta;
                }
            }
        }
    }

    public class ZLevelCoordinates : VerticalDiscretization
    {
        public ZLevelCoordinates(int gridPointsX, int gridPointsZ, double dxi, double deta, double estuaryLength, double estuaryDepth, double[,] x)
            : base(gridPointsX, gridPointsZ, dxi, deta, estuaryLength, estuaryDepth, x)
        {
            InitializeGrid();
        }

        public override void InitializeGrid()
        {
            double dz = estuaryDepth / (gridPointsZ - 1); // Uniform z-level spacing
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    z[i, j] = j * dz; // Z-level: fixed horizontal layers
                }
            }

            ComputeMetricTerms();
        }

        public override void UpdateDepth(double newEstuaryDepth)
        {
            estuaryDepth = newEstuaryDepth;
            InitializeGrid();
        }

        private void ComputeMetricTerms()
        {
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    z_xi[i, j] = i > 0 && i < gridPointsX - 1 ? (z[i + 1, j] - z[i - 1, j]) / (2 * dxi) :
                        i == 0 ? (z[1, j] - z[0, j]) / dxi : (z[gridPointsX - 1, j] - z[gridPointsX - 2, j]) / dxi;
                    z_eta[i, j] = j > 0 && j < gridPointsZ - 1 ? (z[i, j + 1] - z[i, j - 1]) / (2 * deta) :
                        j == 0 ? (z[i, 1] - z[i, 0]) / deta : (z[i, gridPointsZ - 1] - z[i, gridPointsZ - 2]) / deta;
                }
            }
        }
    }
}