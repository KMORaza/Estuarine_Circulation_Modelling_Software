using System;
using System.Collections.Generic;

namespace EstuarineCirculationModeling
{
    public class BifurcatedEstuaryModels
    {
        private readonly int numCellsMain; // Number of cells in main channel
        private readonly int numCellsBranch; // Number of cells per branch
        private readonly double estuaryLength; // Total length of estuary (m)
        private readonly double bifurcationPoint; // Position of bifurcation (fraction of estuary length)
        private readonly double branch1DepthScale; // Depth scaling factor for branch 1
        private readonly double branch2DepthScale; // Depth scaling factor for branch 2
        private readonly double branch1TidalScale; // Tidal velocity scaling for branch 1
        private readonly double branch2TidalScale; // Tidal velocity scaling for branch 2

        public BifurcatedEstuaryModels(
            int numCellsMain = 30,
            int numCellsBranch = 20,
            double estuaryLength = 10000.0,
            double bifurcationPoint = 0.6,
            double branch1DepthScale = 0.8,
            double branch2DepthScale = 1.2,
            double branch1TidalScale = 0.7,
            double branch2TidalScale = 1.3)
        {
            this.numCellsMain = numCellsMain;
            this.numCellsBranch = numCellsBranch;
            this.estuaryLength = estuaryLength;
            this.bifurcationPoint = bifurcationPoint;
            this.branch1DepthScale = branch1DepthScale;
            this.branch2DepthScale = branch2DepthScale;
            this.branch1TidalScale = branch1TidalScale;
            this.branch2TidalScale = branch2TidalScale;
        }

        // Initialize bifurcated grid
        public List<AsymmTidalMix.Cell> InitializeBifurcatedGrid()
        {
            List<AsymmTidalMix.Cell> cells = new List<AsymmTidalMix.Cell>();
            Random rand = new Random(42); // Fixed seed for reproducibility
            double mainChannelLength = estuaryLength * bifurcationPoint;
            double branchLength = estuaryLength * (1.0 - bifurcationPoint);
            double avgCellWidthMain = mainChannelLength / numCellsMain;
            double avgCellWidthBranch = branchLength / numCellsBranch;

            // Create main channel cells (from ocean to bifurcation)
            for (int i = 0; i < numCellsMain; i++)
            {
                double x = i * avgCellWidthMain + rand.NextDouble() * avgCellWidthMain * 0.5;
                double depth = 5.0 + 5.0 * Math.Sin(Math.PI * x / estuaryLength);
                double area = avgCellWidthMain * 100.0 * (0.8 + rand.NextDouble() * 0.4);
                cells.Add(new AsymmTidalMix.Cell(x, depth, area, estuaryLength));
            }

            // Create branch 1 cells
            for (int i = 0; i < numCellsBranch; i++)
            {
                double x = mainChannelLength + i * avgCellWidthBranch + rand.NextDouble() * avgCellWidthBranch * 0.5;
                double depth = (5.0 + 5.0 * Math.Sin(Math.PI * x / estuaryLength)) * branch1DepthScale;
                double area = avgCellWidthBranch * 100.0 * (0.8 + rand.NextDouble() * 0.4);
                cells.Add(new AsymmTidalMix.Cell(x, depth, area, estuaryLength));
            }

            // Create branch 2 cells
            for (int i = 0; i < numCellsBranch; i++)
            {
                double x = mainChannelLength + i * avgCellWidthBranch + rand.NextDouble() * avgCellWidthBranch * 0.5;
                double depth = (5.0 + 5.0 * Math.Sin(Math.PI * x / estuaryLength)) * branch2DepthScale;
                double area = avgCellWidthBranch * 100.0 * (0.8 + rand.NextDouble() * 0.4);
                cells.Add(new AsymmTidalMix.Cell(x, depth, area, estuaryLength));
            }

            // Assign neighbors
            // Main channel: linear connectivity
            for (int i = 0; i < numCellsMain; i++)
            {
                if (i > 0) cells[i].Neighbors.Add(i - 1);
                if (i < numCellsMain - 1) cells[i].Neighbors.Add(i + 1);
                // Bifurcation point connects to first cells of both branches
                if (i == numCellsMain - 1)
                {
                    cells[i].Neighbors.Add(numCellsMain); // First cell of branch 1
                    cells[i].Neighbors.Add(numCellsMain + numCellsBranch); // First cell of branch 2
                }
            }

            // Branch 1: linear connectivity
            for (int i = numCellsMain; i < numCellsMain + numCellsBranch; i++)
            {
                if (i > numCellsMain) cells[i].Neighbors.Add(i - 1); // Connect to previous branch cell
                if (i < numCellsMain + numCellsBranch - 1) cells[i].Neighbors.Add(i + 1); // Connect to next branch cell
                if (i == numCellsMain) cells[i].Neighbors.Add(numCellsMain - 1); // Connect to bifurcation point
            }

            // Branch 2: linear connectivity
            for (int i = numCellsMain + numCellsBranch; i < numCellsMain + 2 * numCellsBranch; i++)
            {
                if (i > numCellsMain + numCellsBranch) cells[i].Neighbors.Add(i - 1); // Connect to previous branch cell
                if (i < numCellsMain + 2 * numCellsBranch - 1) cells[i].Neighbors.Add(i + 1); // Connect to next branch cell
                if (i == numCellsMain + numCellsBranch) cells[i].Neighbors.Add(numCellsMain - 1); // Connect to bifurcation point
            }

            return cells;
        }

        // Determine which branch a cell belongs to (or main channel)
        public int GetBranchIndex(AsymmTidalMix.Cell cell, List<AsymmTidalMix.Cell> cells)
        {
            int index = cells.IndexOf(cell);
            if (index < numCellsMain) return 0; // Main channel
            if (index < numCellsMain + numCellsBranch) return 1; // Branch 1
            return 2; // Branch 2
        }

        // Adjust turbulent viscosity and TKE production based on branch
        public void AdjustMixingDynamics(
            AsymmTidalMix.Cell cell,
            int cellIdx,
            int k,
            double depth,
            double sigmaStep,
            double tidalVelocity,
            ref double nuT,
            ref double kProd,
            List<AsymmTidalMix.Cell> cells)
        {
            int branch = GetBranchIndex(cell, cells);
            double mixingModifier = 1.0;
            double tidalScale = 1.0;

            if (branch == 1) // Branch 1
            {
                mixingModifier = branch1DepthScale; // Shallower branch, less mixing
                tidalScale = branch1TidalScale; // Reduced tidal influence
            }
            else if (branch == 2) // Branch 2
            {
                mixingModifier = branch2DepthScale; // Deeper branch, more mixing
                tidalScale = branch2TidalScale; // Stronger tidal influence
            }

            // Adjust turbulent viscosity
            nuT *= mixingModifier;

            // Adjust TKE production based on tidal velocity and branch depth
            double scaledTidalVelocity = tidalVelocity * tidalScale;
            double shearProduction = nuT * cell.Shear[k] * cell.Shear[k] * mixingModifier;
            kProd += shearProduction * (scaledTidalVelocity / Math.Max(0.01, Math.Abs(tidalVelocity)));
            kProd = Math.Max(0.0, Math.Min(1e-3, kProd)); // Cap for stability
        }

        // Modify tidal velocity for boundary conditions at branch ends
        public double AdjustTidalVelocity(double tidalVelocity, int cellIdx, List<AsymmTidalMix.Cell> cells)
        {
            int branch = GetBranchIndex(cells[cellIdx], cells);
            if (branch == 1)
                return tidalVelocity * branch1TidalScale;
            if (branch == 2)
                return tidalVelocity * branch2TidalScale;
            return tidalVelocity; // Main channel
        }
    }
}