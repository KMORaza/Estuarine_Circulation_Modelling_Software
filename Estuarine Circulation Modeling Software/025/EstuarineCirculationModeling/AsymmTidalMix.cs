using System;
using System.Windows.Forms;
using System.Drawing;
using System.Collections.Generic;

namespace EstuarineCirculationModeling
{
    public class AsymmTidalMix
    {
        private Form atmWindow;
        private Panel visualizationPanel;
        private Panel hodographPanel;
        private TextBox bedShearStressTextBox;
        private TextBox tidalPeriodTextBox;
        private Button startButton;
        private Button pauseButton;
        private Button resetButton;
        private TextBox outputConsoleTextBox;
        private CheckBox showVectorsCheckBox;
        private CheckBox showStreamlinesCheckBox;
        private CheckBox showSalinityCheckBox;
        private CheckBox showDensityCheckBox;
        private CheckBox showTurbidityCheckBox;
        private CheckBox showShearCheckBox;
        private Timer simulationTimer;
        private bool isSimulationRunning;
        private List<Cell> cells; // Unstructured grid cells
        private double bedShearStress; // N/m²
        private double tidalPeriod; // seconds
        private double currentTime;
        private readonly double estuaryLength = 10000.0; // meters
        private readonly int numCells = 50; // Number of horizontal cells
        private readonly static int numSigmaLayers = 10; // Number of vertical sigma layers
        private readonly double dt = 100.0; // seconds, time step
        private readonly double kinematicViscosity = 1e-6; // m²/s
        private readonly double gravity = 9.81; // m/s²
        private readonly double referenceDensity = 1000.0; // kg/m³
        private readonly double coriolisParameter = 1e-4; // s^-1 (mid-latitude approximation)
        private readonly double tidalAmplitude = 1.0; // m/s
        private double tidalPhase;
        private List<(double time, double[] uVelocity, double[] wVelocity)> hodographData; // Store velocity profiles for hodograph

        // Cell structure for unstructured grid
        private class Cell
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

            public Cell(double x, double depth, double area, double estuaryLength)
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

        public AsymmTidalMix()
        {
            InitializeWindow();
            InitializeGrid();
            currentTime = 0.0;
            tidalPhase = 0.0;
            isSimulationRunning = false;
            bedShearStress = 0.1; // Default value
            tidalPeriod = 43200.0; // Default tidal period (12 hours)
            hodographData = new List<(double, double[], double[])>();
            InitializeSimulationTimer();
            ResetSimulation();
        }

        private void InitializeWindow()
        {
            atmWindow = new Form
            {
                Text = "Asymmetric Tidal Mixing (3D Hydrodynamic Solver)",
                Size = new Size(1000, 600),
                FormBorderStyle = FormBorderStyle.FixedDialog,
                MaximizeBox = false,
                Font = new Font("Verdana", 9F)
            };

            // Control Panel
            Panel controlPanel = new Panel
            {
                Location = new Point(10, 10),
                Size = new Size(200, 540),
                BorderStyle = BorderStyle.FixedSingle,
                AutoScroll = true,
                Font = new Font("Verdana", 9F)
            };

            // Bed Shear Stress Input
            Label bedShearStressLabel = new Label
            {
                Location = new Point(10, 10),
                Size = new Size(150, 20),
                Text = "Bed Shear Stress (N/m²):",
                Font = new Font("Verdana", 9F)
            };
            bedShearStressTextBox = new TextBox
            {
                Location = new Point(10, 30),
                Size = new Size(150, 20),
                Text = bedShearStress.ToString("F2"),
                Font = new Font("Verdana", 9F)
            };

            // Tidal Period Input
            Label tidalPeriodLabel = new Label
            {
                Location = new Point(10, 60),
                Size = new Size(150, 20),
                Text = "Tidal Period (s):",
                Font = new Font("Verdana", 9F)
            };
            tidalPeriodTextBox = new TextBox
            {
                Location = new Point(10, 80),
                Size = new Size(150, 20),
                Text = tidalPeriod.ToString("F0"),
                Font = new Font("Verdana", 9F)
            };

            // Visualization Options
            showVectorsCheckBox = new CheckBox
            {
                Location = new Point(10, 110),
                Size = new Size(150, 20),
                Text = "Show Velocity Vectors",
                Checked = true,
                Font = new Font("Verdana", 9F)
            };
            showStreamlinesCheckBox = new CheckBox
            {
                Location = new Point(10, 130),
                Size = new Size(150, 20),
                Text = "Show Streamlines",
                Checked = true,
                Font = new Font("Verdana", 9F)
            };
            showSalinityCheckBox = new CheckBox
            {
                Location = new Point(10, 150),
                Size = new Size(150, 20),
                Text = "Show Salinity Isosurface",
                Checked = true,
                Font = new Font("Verdana", 9F)
            };
            showDensityCheckBox = new CheckBox
            {
                Location = new Point(10, 170),
                Size = new Size(150, 20),
                Text = "Show Density Slice",
                Checked = true,
                Font = new Font("Verdana", 9F)
            };
            showTurbidityCheckBox = new CheckBox
            {
                Location = new Point(10, 190),
                Size = new Size(150, 20),
                Text = "Show Turbidity Fronts",
                Checked = true,
                Font = new Font("Verdana", 9F)
            };
            showShearCheckBox = new CheckBox
            {
                Location = new Point(10, 210),
                Size = new Size(150, 20),
                Text = "Show Shear Layer",
                Checked = true,
                Font = new Font("Verdana", 9F)
            };

            // Buttons
            startButton = new Button
            {
                Location = new Point(10, 240),
                Size = new Size(150, 25),
                Text = "Start",
                FlatStyle = FlatStyle.Flat,
                Font = new Font("Verdana", 9F)
            };
            startButton.Click += StartButton_Click;

            pauseButton = new Button
            {
                Location = new Point(10, 270),
                Size = new Size(150, 25),
                Text = "Pause",
                FlatStyle = FlatStyle.Flat,
                Enabled = false,
                Font = new Font("Verdana", 9F)
            };
            pauseButton.Click += PauseButton_Click;

            resetButton = new Button
            {
                Location = new Point(10, 300),
                Size = new Size(150, 25),
                Text = "Reset",
                FlatStyle = FlatStyle.Flat,
                Font = new Font("Verdana", 9F)
            };
            resetButton.Click += ResetButton_Click;

            // Output Console
            outputConsoleTextBox = new TextBox
            {
                Location = new Point(220, 400),
                Size = new Size(750, 150),
                Multiline = true,
                ReadOnly = true,
                ScrollBars = ScrollBars.Vertical,
                Font = new Font("Verdana", 9F),
                BorderStyle = BorderStyle.FixedSingle
            };

            // Visualization Panel
            visualizationPanel = new Panel
            {
                Location = new Point(220, 10),
                Size = new Size(550, 380),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            visualizationPanel.Paint += VisualizationPanel_Paint;

            // Hodograph Panel
            hodographPanel = new Panel
            {
                Location = new Point(780, 10),
                Size = new Size(190, 380),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            hodographPanel.Paint += HodographPanel_Paint;

            // Add controls to panel
            controlPanel.Controls.AddRange(new Control[] { bedShearStressLabel, bedShearStressTextBox, tidalPeriodLabel, tidalPeriodTextBox,
                showVectorsCheckBox, showStreamlinesCheckBox, showSalinityCheckBox, showDensityCheckBox, showTurbidityCheckBox, showShearCheckBox,
                startButton, pauseButton, resetButton });
            atmWindow.Controls.AddRange(new Control[] { controlPanel, visualizationPanel, hodographPanel, outputConsoleTextBox });

            // Form closing
            atmWindow.FormClosing += (s, e) => simulationTimer?.Stop();
        }

        private void InitializeSimulationTimer()
        {
            simulationTimer = new Timer
            {
                Interval = 100 // Update every 100ms
            };
            simulationTimer.Tick += (s, e) => UpdateSimulation();
        }

        private void StartButton_Click(object sender, EventArgs e)
        {
            try
            {
                bedShearStress = double.Parse(bedShearStressTextBox.Text);
                tidalPeriod = double.Parse(tidalPeriodTextBox.Text);

                // Clamp inputs
                bedShearStress = Math.Max(0.01, Math.Min(1.0, bedShearStress));
                tidalPeriod = Math.Max(3600.0, Math.Min(86400.0, tidalPeriod));

                outputConsoleTextBox.AppendText($"ATM Simulation started with Bed Shear Stress: {bedShearStress:F2} N/m², Tidal Period: {tidalPeriod:F0}s\r\n");
                simulationTimer.Start();
                isSimulationRunning = true;
                UpdateButtonStates();
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Error: {ex.Message}", "Input Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private void PauseButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isSimulationRunning = false;
            outputConsoleTextBox.AppendText("ATM Simulation paused.\r\n");
            UpdateButtonStates();
        }

        private void ResetButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isSimulationRunning = false;
            ResetSimulation();
            outputConsoleTextBox.Clear();
            outputConsoleTextBox.AppendText("ATM Simulation reset.\r\n");
            visualizationPanel.Invalidate();
            hodographPanel.Invalidate();
            UpdateButtonStates();
        }

        private void UpdateButtonStates()
        {
            startButton.Enabled = !isSimulationRunning;
            pauseButton.Enabled = isSimulationRunning;
            resetButton.Enabled = true;
        }

        private void InitializeGrid()
        {
            cells = new List<Cell>();
            Random rand = new Random(42); // Fixed seed for reproducibility

            // Create unstructured grid with variable depth
            double avgCellWidth = estuaryLength / numCells;
            for (int i = 0; i < numCells; i++)
            {
                double x = i * avgCellWidth + rand.NextDouble() * avgCellWidth * 0.5; // Perturb x for unstructured effect
                double depth = 5.0 + 5.0 * Math.Sin(Math.PI * x / estuaryLength); // Variable bathymetry
                double area = avgCellWidth * 100.0 * (0.8 + rand.NextDouble() * 0.4); // Variable cell area
                cells.Add(new Cell(x, depth, area, estuaryLength)); // Pass estuaryLength to constructor
            }

            // Assign neighbors (simplified: each cell connects to previous and next)
            for (int i = 0; i < numCells; i++)
            {
                if (i > 0) cells[i].Neighbors.Add(i - 1);
                if (i < numCells - 1) cells[i].Neighbors.Add(i + 1);
            }
        }

        private void ResetSimulation()
        {
            currentTime = 0.0;
            tidalPhase = 0.0;
            hodographData.Clear();
            foreach (Cell cell in cells)
            {
                for (int k = 0; k < numSigmaLayers; k++)
                {
                    cell.K[k] = 1e-4;
                    cell.Epsilon[k] = 1e-6;
                    cell.U[k] = 0.0;
                    cell.V[k] = 0.0;
                    cell.W[k] = 0.0;
                    cell.Pressure[k] = 0.0;
                    cell.Turbidity[k] = 0.0;
                    cell.Shear[k] = 0.0;
                    double sigma = k / (double)(numSigmaLayers - 1);
                    cell.Salinity[k] = 35.0 * (1.0 - cell.X / estuaryLength); // Linear salinity gradient
                }
            }
        }

        private void UpdateSimulation()
        {
            currentTime += dt;
            tidalPhase = 2 * Math.PI * currentTime / tidalPeriod;

            // Reset hodograph data at the start of each tidal cycle
            if (currentTime >= tidalPeriod)
            {
                hodographData.Clear();
                currentTime = currentTime % tidalPeriod;
                tidalPhase = 2 * Math.PI * currentTime / tidalPeriod;
            }

            // Tidal velocity at boundary (ocean end, i=0)
            double tidalVelocity = tidalAmplitude * Math.Sin(tidalPhase); // Realistic tidal amplitude
            double asymmetryFactor = tidalVelocity >= 0 ? 1.2 : 0.8; // Flood vs. ebb tide

            // Step 1: Compute intermediate velocities (Navier-Stokes with Boussinesq)
            double[,] uStar = new double[numCells, numSigmaLayers];
            double[,] vStar = new double[numCells, numSigmaLayers];
            double[,] wStar = new double[numCells, numSigmaLayers];
            double avgCellWidth = estuaryLength / numCells;

            foreach (Cell cell in cells)
            {
                int i = cells.IndexOf(cell);
                double depth = cell.Depth;
                double sigmaStep = 1.0 / numSigmaLayers;

                for (int k = 0; k < numSigmaLayers; k++)
                {
                    double sigma = k * sigmaStep + sigmaStep / 2;
                    double z = sigma * depth;

                    // Boundary conditions
                    if (i == 0) // Ocean boundary
                    {
                        uStar[i, k] = tidalVelocity * (1.0 - sigma); // Tidal forcing
                        vStar[i, k] = 0.0;
                        wStar[i, k] = 0.0;
                        cell.Salinity[k] = 35.0; // Ocean salinity
                        continue;
                    }
                    else if (i == numCells - 1) // River boundary
                    {
                        uStar[i, k] = 0.1 * (1.0 - sigma); // River inflow
                        vStar[i, k] = 0.0;
                        wStar[i, k] = 0.0;
                        cell.Salinity[k] = 0.0; // Freshwater
                        continue;
                    }

                    // Turbulent viscosity from k-ε model
                    double cMu = 0.09;
                    double nuT = cMu * cell.K[k] * cell.K[k] / cell.Epsilon[k];
                    nuT = Math.Max(1e-6, Math.Min(1e-2, nuT)); // Cap turbulent viscosity

                    // Advection terms (upwind) with stabilization
                    double uAdv = 0.0, vAdv = 0.0, wAdv = 0.0;
                    if (i > 0 && cell.U[k] > 0)
                    {
                        uAdv = cell.U[k] * (cell.U[k] - cells[i - 1].U[k]) / avgCellWidth;
                        vAdv = cell.U[k] * (cell.V[k] - cells[i - 1].V[k]) / avgCellWidth;
                    }
                    else if (i < numCells - 1 && cell.U[k] < 0)
                    {
                        uAdv = cell.U[k] * (cells[i + 1].U[k] - cell.U[k]) / avgCellWidth;
                        vAdv = cell.U[k] * (cells[i + 1].V[k] - cell.V[k]) / avgCellWidth;
                    }
                    if (k < numSigmaLayers - 1 && cell.W[k] > 0)
                    {
                        wAdv = cell.W[k] * (cell.U[k] - cell.U[k + 1]) / (sigmaStep * depth);
                    }
                    else if (k > 0 && cell.W[k] < 0)
                    {
                        wAdv = cell.W[k] * (cell.U[k - 1] - cell.U[k]) / (sigmaStep * depth);
                    }
                    // Cap advection terms
                    uAdv = Math.Max(-10.0, Math.Min(10.0, uAdv));
                    vAdv = Math.Max(-10.0, Math.Min(10.0, vAdv));
                    wAdv = Math.Max(-10.0, Math.Min(10.0, wAdv));

                    // Diffusion terms
                    double uDiff = 0.0, vDiff = 0.0, wDiff = 0.0;
                    foreach (int neighborIdx in cell.Neighbors)
                    {
                        Cell neighbor = cells[neighborIdx];
                        double diff = (kinematicViscosity + nuT) * (neighbor.U[k] - cell.U[k]) / avgCellWidth;
                        uDiff += diff * Math.Min(cell.Area, neighbor.Area) / cell.Volume;
                        diff = (kinematicViscosity + nuT) * (neighbor.V[k] - cell.V[k]) / avgCellWidth;
                        vDiff += diff * Math.Min(cell.Area, neighbor.Area) / cell.Volume;
                    }
                    if (k < numSigmaLayers - 1)
                    {
                        wDiff = (kinematicViscosity + nuT) * (cell.W[k + 1] - cell.W[k]) / (sigmaStep * depth);
                    }
                    // Cap diffusion terms
                    uDiff = Math.Max(-10.0, Math.Min(10.0, uDiff));
                    vDiff = Math.Max(-10.0, Math.Min(10.0, vDiff));
                    wDiff = Math.Max(-10.0, Math.Min(10.0, wDiff));

                    // Coriolis and buoyancy
                    double coriolisU = coriolisParameter * cell.V[k];
                    double coriolisV = -coriolisParameter * cell.U[k];
                    double density = referenceDensity + 0.8 * cell.Salinity[k];
                    double buoyancy = -gravity * (density - referenceDensity) / referenceDensity;

                    // Bed friction
                    double bedFrictionU = (k == 0) ? -bedShearStress * asymmetryFactor * cell.U[k] / referenceDensity : 0.0;
                    double bedFrictionV = (k == 0) ? -bedShearStress * asymmetryFactor * cell.V[k] / referenceDensity : 0.0;

                    // Intermediate velocities with bounds
                    uStar[i, k] = cell.U[k] + dt * (-uAdv + uDiff + coriolisU + bedFrictionU);
                    vStar[i, k] = cell.V[k] + dt * (-vAdv + vDiff + coriolisV + bedFrictionV);
                    wStar[i, k] = cell.W[k] + dt * (-wAdv + wDiff + buoyancy);
                    uStar[i, k] = Math.Max(-2.0, Math.Min(2.0, uStar[i, k]));
                    vStar[i, k] = Math.Max(-2.0, Math.Min(2.0, vStar[i, k]));
                    wStar[i, k] = Math.Max(-2.0, Math.Min(2.0, wStar[i, k]));
                }
            }

            // Step 2: Solve pressure Poisson equation
            double[,] pressure = new double[numCells, numSigmaLayers];
            for (int iter = 0; iter < 20; iter++) // SOR iterations
            {
                for (int i = 0; i < numCells; i++)
                {
                    Cell cell = cells[i];
                    double sigmaStep = 1.0 / numSigmaLayers;
                    for (int k = 0; k < numSigmaLayers; k++)
                    {
                        if (i == 0 || i == numCells - 1) continue; // Skip boundaries
                        double div = 0.0;
                        if (i > 0)
                            div += (uStar[i, k] - uStar[i - 1, k]) / avgCellWidth;
                        if (i < numCells - 1)
                            div += (uStar[i + 1, k] - uStar[i, k]) / avgCellWidth;
                        if (k < numSigmaLayers - 1)
                            div += (wStar[i, k + 1] - wStar[i, k]) / (sigmaStep * cell.Depth);
                        if (k > 0)
                            div += (wStar[i, k] - wStar[i, k - 1]) / (sigmaStep * cell.Depth);
                        div = Math.Max(-10.0, Math.Min(10.0, div)); // Cap divergence

                        double p = 0.0;
                        int count = 0;
                        if (i > 0) { p += pressure[i - 1, k]; count++; }
                        if (i < numCells - 1) { p += pressure[i + 1, k]; count++; }
                        if (k > 0) { p += pressure[i, k - 1]; count++; }
                        if (k < numSigmaLayers - 1) { p += pressure[i, k + 1]; count++; }
                        pressure[i, k] = (1.0 - 1.5) * cell.Pressure[k] + 1.5 * (p / count - div * referenceDensity / dt);
                        pressure[i, k] = Math.Max(-1e5, Math.Min(1e5, pressure[i, k])); // Cap pressure
                    }
                }
            }

            // Step 3: Correct velocities
            foreach (Cell cell in cells)
            {
                int i = cells.IndexOf(cell);
                double sigmaStep = 1.0 / numSigmaLayers;
                for (int k = 0; k < numSigmaLayers; k++)
                {
                    if (i == 0 || i == numCells - 1) continue; // Skip boundaries
                    double dpdx = (i < numCells - 1) ? (pressure[i + 1, k] - pressure[i, k]) / avgCellWidth : 0.0;
                    double dpdz = (k < numSigmaLayers - 1) ? (pressure[i, k + 1] - pressure[i, k]) / (sigmaStep * cell.Depth) : 0.0;
                    dpdx = Math.Max(-1e3, Math.Min(1e3, dpdx)); // Cap pressure gradient
                    dpdz = Math.Max(-1e3, Math.Min(1e3, dpdz));
                    cell.U[k] = uStar[i, k] - dt / referenceDensity * dpdx;
                    cell.V[k] = vStar[i, k];
                    cell.W[k] = wStar[i, k] - dt / referenceDensity * dpdz;
                    // Cap final velocities
                    cell.U[k] = Math.Max(-2.0, Math.Min(2.0, cell.U[k]));
                    cell.V[k] = Math.Max(-2.0, Math.Min(2.0, cell.V[k]));
                    cell.W[k] = Math.Max(-2.0, Math.Min(2.0, cell.W[k]));
                    cell.Pressure[k] = pressure[i, k];

                    // No-slip at bed
                    if (k == 0)
                    {
                        cell.U[k] = 0.0;
                        cell.V[k] = 0.0;
                        cell.W[k] = 0.0;
                    }
                }
            }

            // Step 4: Update turbulence (k-ε model)
            foreach (Cell cell in cells)
            {
                int i = cells.IndexOf(cell);
                double depth = cell.Depth;
                double sigmaStep = 1.0 / numSigmaLayers;

                for (int k = 0; k < numSigmaLayers; k++)
                {
                    double sigma = k * sigmaStep + sigmaStep / 2;
                    double nuT = 0.09 * cell.K[k] * cell.K[k] / cell.Epsilon[k];
                    nuT = Math.Max(1e-6, Math.Min(1e-2, nuT)); // Cap turbulent viscosity

                    // Turbulent kinetic energy (k) equation
                    double kAdv = 0.0, kDiff = 0.0, kProd = 0.0;
                    if (i > 0 && cell.U[k] > 0)
                        kAdv = cell.U[k] * (cell.K[k] - cells[i - 1].K[k]) / avgCellWidth;
                    else if (i < numCells - 1 && cell.U[k] < 0)
                        kAdv = cell.U[k] * (cells[i + 1].K[k] - cell.K[k]) / avgCellWidth;
                    kAdv = Math.Max(-1e-3, Math.Min(1e-3, kAdv));
                    foreach (int neighborIdx in cell.Neighbors)
                    {
                        Cell neighbor = cells[neighborIdx];
                        kDiff += (kinematicViscosity + nuT / 1.0) * (neighbor.K[k] - cell.K[k]) / avgCellWidth * Math.Min(cell.Area, neighbor.Area) / cell.Volume;
                    }
                    kDiff = Math.Max(-1e-3, Math.Min(1e-3, kDiff));
                    kProd = nuT * cell.Shear[k] * cell.Shear[k]; // Shear production
                    kProd = Math.Max(0.0, Math.Min(1e-3, kProd)); // Tighter cap on production
                    cell.K[k] += dt * (-kAdv + kDiff + kProd - cell.Epsilon[k]);
                    cell.K[k] = Math.Max(1e-6, Math.Min(1e-3, cell.K[k])); // Tighter cap on K

                    // Dissipation rate (ε) equation
                    double eAdv = 0.0, eDiff = 0.0;
                    if (i > 0 && cell.U[k] > 0)
                        eAdv = cell.U[k] * (cell.Epsilon[k] - cells[i - 1].Epsilon[k]) / avgCellWidth;
                    else if (i < numCells - 1 && cell.U[k] < 0)
                        eAdv = cell.U[k] * (cells[i + 1].Epsilon[k] - cell.Epsilon[k]) / avgCellWidth;
                    eAdv = Math.Max(-1e-5, Math.Min(1e-5, eAdv));
                    foreach (int neighborIdx in cell.Neighbors)
                    {
                        Cell neighbor = cells[neighborIdx];
                        eDiff += (kinematicViscosity + nuT / 1.3) * (neighbor.Epsilon[k] - cell.Epsilon[k]) / avgCellWidth * Math.Min(cell.Area, neighbor.Area) / cell.Volume;
                    }
                    eDiff = Math.Max(-1e-5, Math.Min(1e-5, eDiff));
                    double eProd = cell.K[k] > 1e-6 ? 1.44 * cell.Epsilon[k] / cell.K[k] * kProd : 0.0;
                    double eDiss = cell.K[k] > 1e-6 ? 1.92 * cell.Epsilon[k] * cell.Epsilon[k] / cell.K[k] : 0.0;
                    cell.Epsilon[k] += dt * (-eAdv + eDiff + eProd - eDiss);
                    cell.Epsilon[k] = Math.Max(1e-8, Math.Min(1e-5, cell.Epsilon[k])); // Tighter cap on ε

                    // Salinity advection and diffusion
                    double salinityChange = 0.0;
                    if (i > 0 && cell.U[k] > 0)
                        salinityChange -= cell.U[k] * (cell.Salinity[k] - cells[i - 1].Salinity[k]) / avgCellWidth;
                    else if (i < numCells - 1 && cell.U[k] < 0)
                        salinityChange -= cell.U[k] * (cells[i + 1].Salinity[k] - cell.Salinity[k]) / avgCellWidth;
                    if (k < numSigmaLayers - 1)
                    {
                        double mixingCoeff = kinematicViscosity + nuT / 1.0;
                        double flux = mixingCoeff * (cell.Salinity[k] - cell.Salinity[k + 1]) / (sigmaStep * depth);
                        salinityChange -= flux * dt / (sigmaStep * depth);
                        if (k > 0)
                        {
                            flux = mixingCoeff * (cell.Salinity[k - 1] - cell.Salinity[k]) / (sigmaStep * depth);
                            salinityChange += flux * dt / (sigmaStep * depth);
                        }
                    }
                    salinityChange = Math.Max(-10.0, Math.Min(10.0, salinityChange)); // Cap salinity change
                    cell.Salinity[k] = Math.Max(0.0, Math.Min(35.0, cell.Salinity[k] + salinityChange)); // Cap salinity

                    // Turbidity (tied to turbulence)
                    double turbidityChange = 0.0;
                    if (i > 0 && cell.U[k] > 0)
                        turbidityChange -= cell.U[k] * (cell.Turbidity[k] - cells[i - 1].Turbidity[k]) / avgCellWidth;
                    else if (i < numCells - 1 && cell.U[k] < 0)
                        turbidityChange -= cell.U[k] * (cells[i + 1].Turbidity[k] - cell.Turbidity[k]) / avgCellWidth;
                    if (k < numSigmaLayers - 1)
                    {
                        double mixingCoeff = kinematicViscosity + nuT / 1.0;
                        double flux = mixingCoeff * (cell.Turbidity[k] - cell.Turbidity[k + 1]) / (sigmaStep * depth);
                        turbidityChange -= flux * dt / (sigmaStep * depth);
                        if (k > 0)
                        {
                            flux = mixingCoeff * (cell.Turbidity[k - 1] - cell.Turbidity[k]) / (sigmaStep * depth);
                            turbidityChange += flux * dt / (sigmaStep * depth);
                        }
                    }
                    turbidityChange = Math.Max(-1.0, Math.Min(1.0, turbidityChange)); // Cap turbidity change
                    double turbidityProduction = (k == 0) ? 0.5 * cell.K[k] : 0.0;
                    cell.Turbidity[k] = Math.Max(0.0, Math.Min(100.0, cell.Turbidity[k] + dt * (turbidityProduction - cell.Turbidity[k] / 1000.0 + turbidityChange))); // Cap turbidity

                    // Shear (magnitude of velocity gradient) with stabilization
                    if (k < numSigmaLayers - 1)
                    {
                        double dz = sigmaStep * depth;
                        double du_dz = (cell.U[k] - cell.U[k + 1]) / dz;
                        double dv_dz = (cell.V[k] - cell.V[k + 1]) / dz;
                        du_dz = Math.Max(-100.0, Math.Min(100.0, du_dz));
                        dv_dz = Math.Max(-100.0, Math.Min(100.0, dv_dz));
                        cell.Shear[k] = Math.Sqrt(du_dz * du_dz + dv_dz * dv_dz);
                        cell.Shear[k] = Math.Min(cell.Shear[k], 10.0);
                    }
                    else
                    {
                        cell.Shear[k] = 0.0;
                    }
                }
            }

            // Store velocity profile for hodograph (mid-estuary cell)
            if ((currentTime % (tidalPeriod / 12)) < dt)
            {
                Cell midCell = cells[numCells / 2];
                double[] uProfile = new double[numSigmaLayers];
                double[] wProfile = new double[numSigmaLayers];
                Array.Copy(midCell.U, uProfile, numSigmaLayers);
                Array.Copy(midCell.W, wProfile, numSigmaLayers);
                hodographData.Add((currentTime, uProfile, wProfile));
            }

            visualizationPanel.Invalidate();
            hodographPanel.Invalidate();
            UpdateOutputConsole();
        }

        private void UpdateOutputConsole()
        {
            double avgTurbulence = 0.0, maxVelocity = 0.0, avgSalinity = 0.0, avgTurbidity = 0.0, maxShear = 0.0, avgK = 0.0;
            int totalLayers = 0;

            foreach (Cell cell in cells)
            {
                for (int k = 0; k < numSigmaLayers; k++)
                {
                    avgTurbulence += cell.K[k];
                    double velocityMag = Math.Sqrt(cell.U[k] * cell.U[k] + cell.V[k] * cell.V[k] + cell.W[k] * cell.W[k]);
                    if (!double.IsNaN(velocityMag) && !double.IsInfinity(velocityMag))
                        maxVelocity = Math.Max(maxVelocity, velocityMag);
                    if (!double.IsNaN(cell.Salinity[k]) && !double.IsInfinity(cell.Salinity[k]))
                        avgSalinity += cell.Salinity[k];
                    if (!double.IsNaN(cell.Turbidity[k]) && !double.IsInfinity(cell.Turbidity[k]))
                        avgTurbidity += cell.Turbidity[k];
                    if (!double.IsNaN(cell.Shear[k]) && !double.IsInfinity(cell.Shear[k]))
                        maxShear = Math.Max(maxShear, cell.Shear[k]);
                    if (!double.IsNaN(cell.K[k]) && !double.IsInfinity(cell.K[k]))
                        avgK += cell.K[k];
                    totalLayers++;
                }
            }
            avgTurbulence = totalLayers > 0 ? avgTurbulence / totalLayers : 0.0;
            avgSalinity = totalLayers > 0 ? avgSalinity / totalLayers : 0.0;
            avgTurbidity = totalLayers > 0 ? avgTurbidity / totalLayers : 0.0;
            avgK = totalLayers > 0 ? avgK / totalLayers : 0.0;

            // Log warnings for invalid values
            if (double.IsNaN(maxVelocity) || double.IsInfinity(maxVelocity))
                outputConsoleTextBox.AppendText($"Warning: Invalid Max Velocity at time {currentTime:F2}s\r\n");
            if (double.IsNaN(avgSalinity) || double.IsInfinity(avgSalinity))
                outputConsoleTextBox.AppendText($"Warning: Invalid Avg Salinity at time {currentTime:F2}s\r\n");
            if (double.IsNaN(avgTurbidity) || double.IsInfinity(avgTurbidity))
                outputConsoleTextBox.AppendText($"Warning: Invalid Avg Turbidity at time {currentTime:F2}s\r\n");
            if (double.IsNaN(maxShear) || double.IsInfinity(maxShear))
                outputConsoleTextBox.AppendText($"Warning: Invalid Max Shear at time {currentTime:F2}s\r\n");

            outputConsoleTextBox.AppendText($"Time: {currentTime:F2}s | Avg TKE: {avgK:F6} m²/s² | Max Velocity: {maxVelocity:F4} m/s | Avg Salinity: {avgSalinity:F2} PSU | Avg Turbidity: {avgTurbidity:F4} | Max Shear: {maxShear:F4} 1/s | Tide: {(tidalPhase >= 0 && tidalPhase < Math.PI ? "Flood" : "Ebb")}\r\n");
        }

        private void VisualizationPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            int width = visualizationPanel.Width;
            int height = visualizationPanel.Height;

            // Draw background
            g.Clear(Color.White);

            // Draw tidal phase indicator
            string phaseText = tidalPhase >= 0 && tidalPhase < Math.PI ? "Flood" : "Ebb";
            using (Brush textBrush = new SolidBrush(tidalPhase >= 0 && tidalPhase < Math.PI ? Color.Red : Color.Blue))
            {
                g.DrawString($"Tidal Phase: {phaseText}", new Font("Verdana", 10), textBrush, 10, 10);
            }

            // Compute max values for scaling
            double maxTurbulence = 0.0, maxVelocity = 0.0, maxSalinity = 0.0, maxTurbidity = 0.0, maxShear = 0.0;
            foreach (Cell cell in cells)
            {
                for (int k = 0; k < numSigmaLayers; k++)
                {
                    if (!double.IsNaN(cell.K[k]) && !double.IsInfinity(cell.K[k]))
                        maxTurbulence = Math.Max(maxTurbulence, cell.K[k]);
                    double velocityMag = Math.Sqrt(cell.U[k] * cell.U[k] + cell.V[k] * cell.V[k] + cell.W[k] * cell.W[k]);
                    if (!double.IsNaN(velocityMag) && !double.IsInfinity(velocityMag))
                        maxVelocity = Math.Max(maxVelocity, velocityMag);
                    if (!double.IsNaN(cell.Salinity[k]) && !double.IsInfinity(cell.Salinity[k]))
                        maxSalinity = Math.Max(maxSalinity, cell.Salinity[k]);
                    if (!double.IsNaN(cell.Turbidity[k]) && !double.IsInfinity(cell.Turbidity[k]))
                        maxTurbidity = Math.Max(maxTurbidity, cell.Turbidity[k]);
                    if (!double.IsNaN(cell.Shear[k]) && !double.IsInfinity(cell.Shear[k]))
                        maxShear = Math.Max(maxShear, cell.Shear[k]);
                }
            }
            maxTurbulence = Math.Max(maxTurbulence, 1e-6);
            maxTurbulence = Math.Min(maxTurbulence, 1e-3); // Tighter cap
            maxVelocity = Math.Max(maxVelocity, 0.01);
            maxVelocity = Math.Min(maxVelocity, 2.0); // Cap at realistic velocity
            maxSalinity = Math.Max(maxSalinity, 0.01);
            maxSalinity = Math.Min(maxSalinity, 35.0); // Cap at ocean salinity
            maxTurbidity = Math.Max(maxTurbidity, 0.01);
            maxTurbidity = Math.Min(maxTurbidity, 100.0); // Cap turbidity
            maxShear = Math.Max(maxShear, 0.01);
            maxShear = Math.Min(maxShear, 10.0);

            // Draw density slice (x-z plane) if checked
            if (showDensityCheckBox.Checked)
            {
                for (int i = 0; i < numCells; i++)
                {
                    Cell cell = cells[i];
                    float x = (float)(cell.X / estuaryLength * width);
                    for (int k = 0; k < numSigmaLayers; k++)
                    {
                        double sigma = k / (double)(numSigmaLayers - 1);
                        float y = height - (float)(sigma * height);
                        float cellWidth = (float)(estuaryLength / numCells / estuaryLength * width);
                        float cellHeight = (float)(height / numSigmaLayers);
                        double density = referenceDensity + 0.8 * cell.Salinity[k];
                        if (double.IsNaN(density) || double.IsInfinity(density)) continue;
                        int colorValue = (int)((density - referenceDensity) / (0.8 * maxSalinity) * 255);
                        colorValue = Math.Max(0, Math.Min(255, colorValue));
                        if (float.IsNaN(x) || float.IsInfinity(x) || float.IsNaN(y) || float.IsInfinity(y)) continue;
                        x = Math.Max(0, Math.Min(width, x));
                        y = Math.Max(0, Math.Min(height, y));
                        using (Brush brush = new SolidBrush(Color.FromArgb(colorValue, 0, 255 - colorValue)))
                        {
                            g.FillRectangle(brush, x, y - cellHeight, cellWidth, cellHeight);
                        }
                    }
                }
            }

            // Draw salinity isosurface (15 PSU) with time-lapse effect
            if (showSalinityCheckBox.Checked)
            {
                double isoValue = 15.0; // Isosurface at 15 PSU
                for (int i = 0; i < numCells - 1; i++)
                {
                    Cell cell1 = cells[i];
                    Cell cell2 = cells[i + 1];
                    float x1 = (float)(cell1.X / estuaryLength * width);
                    float x2 = (float)(cell2.X / estuaryLength * width);
                    if (float.IsNaN(x1) || float.IsInfinity(x1) || float.IsNaN(x2) || float.IsInfinity(x2)) continue;
                    x1 = Math.Max(0, Math.Min(width, x1));
                    x2 = Math.Max(0, Math.Min(width, x2));
                    for (int k = 0; k < numSigmaLayers - 1; k++)
                    {
                        double s1 = cell1.Salinity[k];
                        double s2 = cell2.Salinity[k];
                        double s1_next = cell1.Salinity[k + 1];
                        double s2_next = cell2.Salinity[k + 1];
                        if (double.IsNaN(s1) || double.IsInfinity(s1) || double.IsNaN(s2) || double.IsInfinity(s2) ||
                            double.IsNaN(s1_next) || double.IsInfinity(s1_next) || double.IsNaN(s2_next) || double.IsInfinity(s2_next)) continue;
                        double sigma = k / (double)(numSigmaLayers - 1);
                        double sigma_next = (k + 1) / (double)(numSigmaLayers - 1);
                        float y1 = height - (float)(sigma * height);
                        float y2 = height - (float)(sigma_next * height);
                        if (float.IsNaN(y1) || float.IsInfinity(y1) || float.IsNaN(y2) || float.IsInfinity(y2)) continue;
                        y1 = Math.Max(0, Math.Min(height, y1));
                        y2 = Math.Max(0, Math.Min(height, y2));
                        if ((s1 - isoValue) * (s2 - isoValue) <= 0 || (s1 - isoValue) * (s1_next - isoValue) <= 0)
                        {
                            using (Pen pen = new Pen(Color.FromArgb(tidalPhase >= 0 && tidalPhase < Math.PI ? 255 : 128, 0, 255, 0), tidalPhase >= 0 && tidalPhase < Math.PI ? 2f : 1f))
                            {
                                g.DrawLine(pen, x1, y1, x2, y1);
                            }
                        }
                    }
                }
            }

            // Draw shear layer
            if (showShearCheckBox.Checked)
            {
                for (int i = 0; i < numCells - 1; i++)
                {
                    Cell cell1 = cells[i];
                    Cell cell2 = cells[i + 1];
                    float x1 = (float)(cell1.X / estuaryLength * width);
                    float x2 = (float)(cell2.X / estuaryLength * width);
                    if (float.IsNaN(x1) || float.IsInfinity(x1) || float.IsNaN(x2) || float.IsInfinity(x2)) continue;
                    x1 = Math.Max(0, Math.Min(width, x1));
                    x2 = Math.Max(0, Math.Min(width, x2));
                    for (int k = 0; k < numSigmaLayers - 1; k++)
                    {
                        double shear1 = cell1.Shear[k];
                        if (double.IsNaN(shear1) || double.IsInfinity(shear1)) continue;
                        float y = height - (float)(k / (double)(numSigmaLayers - 1) * height);
                        float shearWidth = (float)Math.Min(Math.Max(shear1 / maxShear * 3, 0.5), 10.0);
                        if (float.IsNaN(y) || float.IsInfinity(y) || float.IsNaN(shearWidth) || float.IsInfinity(shearWidth)) continue;
                        y = Math.Max(0, Math.Min(height, y));
                        using (Pen pen = new Pen(Color.FromArgb(tidalPhase >= 0 && tidalPhase < Math.PI ? 255 : 128, 255, 165, 0), shearWidth))
                        {
                            g.DrawLine(pen, x1, y, x2, y);
                        }
                    }
                }
            }

            // Draw turbidity fronts
            if (showTurbidityCheckBox.Checked)
            {
                double threshold = 0.5 * maxTurbidity;
                for (int i = 0; i < numCells - 1; i++)
                {
                    Cell cell1 = cells[i];
                    Cell cell2 = cells[i + 1];
                    float x1 = (float)(cell1.X / estuaryLength * width);
                    float x2 = (float)(cell2.X / estuaryLength * width);
                    if (float.IsNaN(x1) || float.IsInfinity(x1) || float.IsNaN(x2) || float.IsInfinity(x2)) continue;
                    x1 = Math.Max(0, Math.Min(width, x1));
                    x2 = Math.Max(0, Math.Min(width, x2));
                    double avgTurb1 = 0.0, avgTurb2 = 0.0;
                    int validLayers = 0;
                    for (int k = 0; k < numSigmaLayers; k++)
                    {
                        if (!double.IsNaN(cell1.Turbidity[k]) && !double.IsInfinity(cell1.Turbidity[k]))
                            avgTurb1 += cell1.Turbidity[k];
                        if (!double.IsNaN(cell2.Turbidity[k]) && !double.IsInfinity(cell2.Turbidity[k]))
                            avgTurb2 += cell2.Turbidity[k];
                        validLayers++;
                    }
                    if (validLayers == 0) continue;
                    avgTurb1 /= validLayers;
                    avgTurb2 /= validLayers;
                    if (double.IsNaN(avgTurb1) || double.IsInfinity(avgTurb1) || double.IsNaN(avgTurb2) || double.IsInfinity(avgTurb2)) continue;
                    if ((avgTurb1 - threshold) * (avgTurb2 - threshold) <= 0)
                    {
                        float y = height - (float)(0.5 * height);
                        if (float.IsNaN(y) || float.IsInfinity(y)) continue;
                        y = Math.Max(0, Math.Min(height, y));
                        using (Pen pen = new Pen(Color.FromArgb(tidalPhase >= 0 && tidalPhase < Math.PI ? 255 : 128, 139, 69, 19), tidalPhase >= 0 && tidalPhase < Math.PI ? 2f : 1f))
                        {
                            g.DrawLine(pen, x1, y, x2, y);
                        }
                    }
                }
            }

            // Draw turbulence and velocity profiles
            for (int i = 0; i < numCells - 1; i++)
            {
                Cell cell1 = cells[i];
                Cell cell2 = cells[i + 1];
                double avgTurb1 = 0.0, avgTurb2 = 0.0;
                double avgVel1 = 0.0, avgVel2 = 0.0;
                int validLayers = 0;
                for (int k = 0; k < numSigmaLayers; k++)
                {
                    if (!double.IsNaN(cell1.K[k]) && !double.IsInfinity(cell1.K[k]))
                        avgTurb1 += cell1.K[k];
                    if (!double.IsNaN(cell2.K[k]) && !double.IsInfinity(cell2.K[k]))
                        avgTurb2 += cell2.K[k];
                    double vel1 = Math.Sqrt(cell1.U[k] * cell1.U[k] + cell1.V[k] * cell1.V[k]);
                    double vel2 = Math.Sqrt(cell2.U[k] * cell2.U[k] + cell2.V[k] * cell2.V[k]);
                    if (!double.IsNaN(vel1) && !double.IsInfinity(vel1))
                        avgVel1 += vel1;
                    if (!double.IsNaN(vel2) && !double.IsInfinity(vel2))
                        avgVel2 += vel2;
                    validLayers++;
                }
                if (validLayers == 0) continue;
                avgTurb1 /= validLayers;
                avgTurb2 /= validLayers;
                avgVel1 /= validLayers;
                avgVel2 /= validLayers;
                if (double.IsNaN(avgTurb1) || double.IsInfinity(avgTurb1) || double.IsNaN(avgTurb2) || double.IsInfinity(avgTurb2) ||
                    double.IsNaN(avgVel1) || double.IsInfinity(avgVel1) || double.IsNaN(avgVel2) || double.IsInfinity(avgVel2)) continue;
                float x1 = (float)(cell1.X / estuaryLength * width);
                float x2 = (float)(cell2.X / estuaryLength * width);
                float y1Turb = height - (float)(avgTurb1 / maxTurbulence * height / 2);
                float y2Turb = height - (float)(avgTurb2 / maxTurbulence * height / 2);
                float y1Vel = height / 2 - (float)(avgVel1 / maxVelocity * height / 4);
                float y2Vel = height / 2 - (float)(avgVel2 / maxVelocity * height / 4);
                if (float.IsNaN(x1) || float.IsInfinity(x1) || float.IsNaN(x2) || float.IsInfinity(x2) ||
                    float.IsNaN(y1Turb) || float.IsInfinity(y1Turb) || float.IsNaN(y2Turb) || float.IsInfinity(y2Turb) ||
                    float.IsNaN(y1Vel) || float.IsInfinity(y1Vel) || float.IsNaN(y2Vel) || float.IsInfinity(y2Vel)) continue;
                x1 = Math.Max(0, Math.Min(width, x1));
                x2 = Math.Max(0, Math.Min(width, x2));
                y1Turb = Math.Max(0, Math.Min(height, y1Turb));
                y2Turb = Math.Max(0, Math.Min(height, y2Turb));
                y1Vel = Math.Max(0, Math.Min(height, y1Vel));
                y2Vel = Math.Max(0, Math.Min(height, y2Vel));
                using (Pen turbPen = new Pen(Color.FromArgb(tidalPhase >= 0 && tidalPhase < Math.PI ? 255 : 128, 0, 0, 255), tidalPhase >= 0 && tidalPhase < Math.PI ? 2f : 1f))
                {
                    g.DrawLine(turbPen, x1, y1Turb, x2, y2Turb);
                }
                using (Pen velPen = new Pen(Color.FromArgb(tidalPhase >= 0 && tidalPhase < Math.PI ? 255 : 128, 255, 0, 0), tidalPhase >= 0 && tidalPhase < Math.PI ? 2f : 1f))
                {
                    g.DrawLine(velPen, x1, y1Vel, x2, y2Vel);
                }
            }

            // Draw velocity vectors
            if (showVectorsCheckBox.Checked)
            {
                for (int i = 0; i < numCells; i += 5)
                {
                    Cell cell = cells[i];
                    float x = (float)(cell.X / estuaryLength * width);
                    for (int k = 0; k < numSigmaLayers; k += 2)
                    {
                        double sigma = k / (double)(numSigmaLayers - 1);
                        float y = height - (float)(sigma * height);
                        float vx = (float)(cell.U[k] / maxVelocity * 20);
                        float vz = (float)(cell.W[k] / maxVelocity * 20);
                        if (float.IsNaN(x) || float.IsInfinity(x) || float.IsNaN(y) || float.IsInfinity(y) ||
                            float.IsNaN(vx) || float.IsInfinity(vx) || float.IsNaN(vz) || float.IsInfinity(vz)) continue;
                        x = Math.Max(0, Math.Min(width, x));
                        y = Math.Max(0, Math.Min(height, y));
                        using (Pen pen = new Pen(Color.FromArgb(tidalPhase >= 0 && tidalPhase < Math.PI ? 255 : 128, 0, 0, 0), tidalPhase >= 0 && tidalPhase < Math.PI ? 2f : 1f))
                        {
                            g.DrawLine(pen, x, y, x + vx, y + vz);
                            g.DrawLine(pen, x + vx, y + vz, x + vx - 3, y + vz - 3);
                            g.DrawLine(pen, x + vx, y + vz, x + vx - 3, y + vz + 3);
                        }
                    }
                }
            }

            // Draw streamlines
            if (showStreamlinesCheckBox.Checked)
            {
                for (int start = 0; start < numCells; start += 10)
                {
                    float x = (float)(cells[start].X / estuaryLength * width);
                    float y = height / 2;
                    List<PointF> streamline = new List<PointF> { new PointF(x, y) };
                    for (int step = 0; step < 20; step++)
                    {
                        int cellIdx = -1;
                        double minDist = double.MaxValue;
                        for (int i = 0; i < numCells; i++)
                        {
                            double dist = Math.Abs(cells[i].X - (x / width * estuaryLength));
                            if (dist < minDist)
                            {
                                minDist = dist;
                                cellIdx = i;
                            }
                        }
                        if (cellIdx < 0 || cellIdx >= numCells) break;
                        Cell cell = cells[cellIdx];
                        double avgU = 0.0;
                        int validLayers = 0;
                        for (int k = 0; k < numSigmaLayers; k++)
                        {
                            if (!double.IsNaN(cell.U[k]) && !double.IsInfinity(cell.U[k]))
                            {
                                avgU += cell.U[k];
                                validLayers++;
                            }
                        }
                        if (validLayers == 0) break;
                        avgU /= validLayers;
                        if (double.IsNaN(avgU) || double.IsInfinity(avgU)) break;
                        x += (float)(avgU / maxVelocity * 20);
                        if (float.IsNaN(x) || float.IsInfinity(x) || x >= width || x < 0) break;
                        x = Math.Max(0, Math.Min(width, x));
                        streamline.Add(new PointF(x, y));
                    }
                    if (streamline.Count > 1)
                    {
                        using (Pen pen = new Pen(Color.FromArgb(tidalPhase >= 0 && tidalPhase < Math.PI ? 255 : 128, 128, 0, 128), tidalPhase >= 0 && tidalPhase < Math.PI ? 2f : 1f))
                        {
                            g.DrawLines(pen, streamline.ToArray());
                        }
                    }
                }
            }
        }

        private void HodographPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            int width = hodographPanel.Width;
            int height = hodographPanel.Height;

            // Draw background and axes
            g.Clear(Color.White);
            g.DrawLine(Pens.Black, width / 2, 0, width / 2, height);
            g.DrawLine(Pens.Black, 0, height / 2, width, height / 2);

            // Draw hodograph
            double maxVelocity = 0.0;
            foreach (var (time, uVelocity, wVelocity) in hodographData)
            {
                for (int k = 0; k < numSigmaLayers; k++)
                {
                    double velocityMag = Math.Sqrt(uVelocity[k] * uVelocity[k] + wVelocity[k] * wVelocity[k]);
                    if (!double.IsNaN(velocityMag) && !double.IsInfinity(velocityMag))
                        maxVelocity = Math.Max(maxVelocity, velocityMag);
                }
            }
            maxVelocity = Math.Max(maxVelocity, 0.01);
            maxVelocity = Math.Min(maxVelocity, 2.0); // Cap at realistic velocity

            for (int k = 0; k < numSigmaLayers; k += 2)
            {
                List<PointF> hodographPoints = new List<PointF>();
                foreach (var (time, uVelocity, wVelocity) in hodographData)
                {
                    float vx = (float)(uVelocity[k] / maxVelocity * (width / 4));
                    float vz = (float)(wVelocity[k] / maxVelocity * (height / 4));
                    float x = width / 2 + vx;
                    float y = height / 2 + vz;
                    if (float.IsNaN(x) || float.IsInfinity(x) || float.IsNaN(y) || float.IsInfinity(y)) continue;
                    x = Math.Max(0, Math.Min(width, x));
                    y = Math.Max(0, Math.Min(height, y));
                    hodographPoints.Add(new PointF(x, y));
                }
                if (hodographPoints.Count > 1)
                {
                    using (Pen pen = new Pen(Color.FromArgb(tidalPhase >= 0 && tidalPhase < Math.PI ? 255 : 128, 0, 0, 255), 1f))
                    {
                        g.DrawLines(pen, hodographPoints.ToArray());
                    }
                    PointF current = hodographPoints[hodographPoints.Count - 1];
                    g.FillEllipse(Brushes.Red, current.X - 3, current.Y - 3, 6, 6);
                }
            }
        }

        public void ShowAsymmTidalMixWindow()
        {
            if (!atmWindow.Visible)
                atmWindow.Show();
            else
                atmWindow.BringToFront();
        }
    }
}