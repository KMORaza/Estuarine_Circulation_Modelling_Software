using System;
using System.Collections.Generic;

namespace EstuarineCirculationModeling
{
    public class Cell
    {
        public double X; // Representative x-coordinate (m)
        public double Depth; // Local water depth (m)
        public double Volume; // Cell volume (m³)
        public double[] K; // Turbulent kinetic energy per sigma layer (m²/s²)
        public double[] Epsilon; // Turbulence dissipation rate per sigma layer (m²/s³)
        public double[] U; // x-velocity per sigma layer (m/s)
        public double[] V; // y-velocity per sigma layer (m/s)
        public double[] W; // z-velocity per sigma layer (m/s)
        public double[] Pressure; // Pressure per sigma layer (Pa)
        public double[] Salinity; // Salinity per sigma layer (PSU)
        public double[] Turbidity; // Turbidity per sigma layer (arbitrary units)
        public double[] Shear; // Shear (velocity gradient) per sigma layer (1/s)
        public List<int> Neighbors; // Indices of neighboring cells
        public double Area; // Horizontal area of cell (m²)

        public Cell(double x, double depth, double area, double estuaryLength, int numSigmaLayers)
        {
            X = x;
            Depth = depth;
            Area = area;
            Volume = area * depth;
            K = new double[numSigmaLayers];
            Epsilon = new double[numSigmaLayers];
            U = new double[numSigmaLayers];
            V = new double[numSigmaLayers];
            W = new double[numSigmaLayers];
            Pressure = new double[numSigmaLayers];
            Salinity = new double[numSigmaLayers];
            Turbidity = new double[numSigmaLayers];
            Shear = new double[numSigmaLayers];
            Neighbors = new List<int>();
            // Initialize salinity (35 PSU at ocean, 0 PSU at river)
            for (int k = 0; k < numSigmaLayers; k++)
            {
                double sigma = k / (double)(numSigmaLayers - 1);
                Salinity[k] = 35.0 * (1.0 - x / estuaryLength); // Linear salinity gradient
                K[k] = 1e-4; // Initial turbulent kinetic energy
                Epsilon[k] = 1e-6; // Initial dissipation rate
            }
        }
    }
}