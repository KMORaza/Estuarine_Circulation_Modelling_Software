using System;
using System.Drawing;
using System.Windows.Forms;

namespace EstuarineCirculationModeling
{
    public static class LargeEddySim
    {
        private static int gridPointsX = 50; // Grid points in x-direction
        private static int gridPointsZ = 20; // Grid points in z-direction
        private static double estuaryLength = 10000.0; // m
        private static double estuaryDepth = 10.0; // m
        private static double dx; // Spatial step x (m)
        private static double dz; // Spatial step z (m)
        private static double dt = 1.0; // Time step (s)
        private static double smagorinskyConstant = 0.1; // Smagorinsky constant for LES
        private static double kinematicViscosity = 1e-6; // m²/s
        private static double salinityOcean = 35.0; // PSU
        private static double temperatureOcean = 20.0; // °C
        private static double riverInflow = 0.1; // m³/s
        private static double tidalAmplitude = 1.0; // m
        private static double tidalPeriod = 43200.0; // s (12 hours)
        private static double simulationTime = 0.0; // Current simulation time (s)
        private static double[,] velocityU; // Horizontal velocity (m/s)
        private static double[,] velocityW; // Vertical velocity (m/s)
        private static double[,] salinity; // Salinity (PSU)
        private static double[,] temperature; // Temperature (°C)
        private static double[,] vorticity; // Vorticity (s^-1)
        private static double[,] eddyViscosity; // Eddy viscosity (m²/s)
        private static readonly double g = 9.81; // Gravitational acceleration (m/s²)
        private static readonly double densityReference = 1000.0; // Reference density (kg/m³)
        private static Form lesForm;
        private static Panel simulationPanel;
        private static Timer simulationTimer;
        private static Label statusLabel;
        private static bool isSimulationRunning;

        public static void ShowLargeEddyWindow()
        {
            // Initialize simulation parameters
            InitializeSimulation();
            isSimulationRunning = false;

            // Create form
            lesForm = new Form
            {
                Text = "Large Eddy Simulation",
                ClientSize = new Size(360, 640), // Mobile-like size (9:16 aspect ratio)
                FormBorderStyle = FormBorderStyle.FixedSingle,
                MaximizeBox = false,
                StartPosition = FormStartPosition.CenterScreen,
                BackColor = Color.White,
                ForeColor = Color.Black,
                Font = new Font("Segoe UI", 10F)
            };

            // Setup UI
            SetupUI();

            // Update status
            UpdateStatus("Ready");

            // Show form
            lesForm.ShowDialog();
        }

        private static void InitializeSimulation()
        {
            // Ensure valid grid sizes
            if (gridPointsX < 20 || gridPointsX > 100 || gridPointsZ < 10 || gridPointsZ > 50)
            {
                gridPointsX = 50; // Default fallback
                gridPointsZ = 20;
                UpdateStatus("Grid sizes reset to default (50x20)");
            }

            try
            {
                dx = estuaryLength / gridPointsX;
                dz = estuaryDepth / gridPointsZ;
                velocityU = new double[gridPointsX, gridPointsZ];
                velocityW = new double[gridPointsX, gridPointsZ];
                salinity = new double[gridPointsX, gridPointsZ];
                temperature = new double[gridPointsX, gridPointsZ];
                vorticity = new double[gridPointsX, gridPointsZ];
                eddyViscosity = new double[gridPointsX, gridPointsZ];
                ResetSimulation();
            }
            catch (Exception ex)
            {
                UpdateStatus($"Initialization failed: {ex.Message}");
                MessageBox.Show($"Error initializing simulation: {ex.Message}", "Initialization Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private static void SetupUI()
        {
            // Border panel to simulate mobile screen frame
            Panel borderPanel = new Panel
            {
                Location = new Point(5, 5),
                Size = new Size(350, 630),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            lesForm.Controls.Add(borderPanel);

            // Control panel
            Panel controlPanel = new Panel
            {
                Location = new Point(10, 10),
                Size = new Size(330, 260),
                AutoScroll = true, // Scrollable for future controls
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White,
                ForeColor = Color.Black
            };
            borderPanel.Controls.Add(controlPanel);

            // Simulation panel
            simulationPanel = new Panel
            {
                Location = new Point(10, 280),
                Size = new Size(330, 340),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            simulationPanel.Paint += (s, e) => DrawSimulation(e.Graphics, simulationPanel.Width, simulationPanel.Height);
            borderPanel.Controls.Add(simulationPanel);

            // Smagorinsky constant
            Label smagorinskyLabel = new Label
            {
                Text = "Smagorinsky Constant:",
                Location = new Point(10, 10),
                AutoSize = true,
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black
            };
            TextBox smagorinskyTextBox = new TextBox
            {
                Location = new Point(10, 35),
                Size = new Size(310, 25),
                Text = smagorinskyConstant.ToString("F2"),
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            smagorinskyTextBox.TextChanged += (s, e) => ValidateSmagorinsky(smagorinskyTextBox);

            // Grid resolution
            Label gridLabel = new Label
            {
                Text = "Grid Points (XxZ):",
                Location = new Point(10, 70),
                AutoSize = true,
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black
            };
            TextBox gridTextBox = new TextBox
            {
                Location = new Point(10, 95),
                Size = new Size(310, 25),
                Text = $"{gridPointsX}x{gridPointsZ}",
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            gridTextBox.TextChanged += (s, e) => ValidateGrid(gridTextBox);

            // River inflow
            Label inflowLabel = new Label
            {
                Text = "River Inflow (m³/s):",
                Location = new Point(10, 130),
                AutoSize = true,
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black
            };
            TextBox inflowTextBox = new TextBox
            {
                Location = new Point(10, 155),
                Size = new Size(310, 25),
                Text = riverInflow.ToString("F2"),
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            inflowTextBox.TextChanged += (s, e) => ValidateInflow(inflowTextBox);

            // Tidal amplitude
            Label tidalLabel = new Label
            {
                Text = "Tidal Amplitude (m):",
                Location = new Point(10, 190),
                AutoSize = true,
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black
            };
            TextBox tidalTextBox = new TextBox
            {
                Location = new Point(10, 215),
                Size = new Size(310, 25),
                Text = tidalAmplitude.ToString("F2"),
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            tidalTextBox.TextChanged += (s, e) => ValidateTidal(tidalTextBox);

            // Apply button
            Button applyButton = new Button
            {
                Text = "Apply",
                Location = new Point(10, 250),
                Size = new Size(100, 30),
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = Color.FromArgb(200, 200, 200),
                ForeColor = Color.Black
            };
            applyButton.Click += (s, e) =>
            {
                if (ValidateAllInputs(smagorinskyTextBox, gridTextBox, inflowTextBox, tidalTextBox))
                {
                    InitializeSimulation();
                    simulationPanel.Invalidate();
                    UpdateStatus("Parameters applied");
                }
            };

            // Close button
            Button closeButton = new Button
            {
                Text = "Close",
                Location = new Point(120, 250),
                Size = new Size(100, 30),
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = Color.FromArgb(200, 200, 200),
                ForeColor = Color.Black
            };
            closeButton.Click += (s, e) =>
            {
                simulationTimer.Stop();
                lesForm.Close();
            };

            // Start button
            Button startButton = new Button
            {
                Text = "Start",
                Location = new Point(10, 290),
                Size = new Size(100, 30),
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = Color.FromArgb(200, 200, 200),
                ForeColor = Color.Black
            };
            startButton.Click += (s, e) =>
            {
                if (ValidateAllInputs(smagorinskyTextBox, gridTextBox, inflowTextBox, tidalTextBox))
                {
                    if (salinity == null || temperature == null || salinity.GetLength(0) != gridPointsX || salinity.GetLength(1) != gridPointsZ)
                    {
                        InitializeSimulation();
                    }
                    simulationTimer.Start();
                    isSimulationRunning = true;
                    UpdateStatus("Running");
                }
            };

            // Pause button
            Button pauseButton = new Button
            {
                Text = "Pause",
                Location = new Point(120, 290),
                Size = new Size(100, 30),
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = Color.FromArgb(200, 200, 200),
                ForeColor = Color.Black
            };
            pauseButton.Click += (s, e) =>
            {
                simulationTimer.Stop();
                isSimulationRunning = false;
                UpdateStatus("Paused");
            };

            // Reset button
            Button resetButton = new Button
            {
                Text = "Reset",
                Location = new Point(230, 290),
                Size = new Size(100, 30),
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = Color.FromArgb(200, 200, 200),
                ForeColor = Color.Black
            };
            resetButton.Click += (s, e) =>
            {
                simulationTimer.Stop();
                isSimulationRunning = false;
                InitializeSimulation();
                simulationPanel.Invalidate();
                UpdateStatus("Reset");
            };

            // Status label
            statusLabel = new Label
            {
                Text = "Ready",
                Location = new Point(10, 330),
                Size = new Size(310, 40),
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black
            };

            // Simulation timer
            simulationTimer = new Timer
            {
                Interval = 100 // Update every 100ms
            };
            simulationTimer.Tick += (s, e) =>
            {
                if (salinity == null || temperature == null || salinity.GetLength(0) != gridPointsX || salinity.GetLength(1) != gridPointsZ)
                {
                    simulationTimer.Stop();
                    isSimulationRunning = false;
                    UpdateStatus("Error: Arrays not initialized");
                    MessageBox.Show("Simulation arrays not properly initialized", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                    return;
                }
                UpdateSimulation();
                simulationPanel.Invalidate();
                UpdateStatus($"Running (Time: {simulationTime:F1} s)");
            };

            // Add controls to control panel
            controlPanel.Controls.AddRange(new Control[] { smagorinskyLabel, smagorinskyTextBox, gridLabel, gridTextBox, inflowLabel, inflowTextBox, tidalLabel, tidalTextBox, applyButton, closeButton, startButton, pauseButton, resetButton, statusLabel });

            // Add panels to form
            lesForm.Controls.Add(borderPanel);
        }

        private static bool ValidateSmagorinsky(TextBox textBox)
        {
            if (double.TryParse(textBox.Text, out double value))
            {
                smagorinskyConstant = Math.Max(0.01, Math.Min(0.5, value));
                textBox.Text = smagorinskyConstant.ToString("F2");
                UpdateStatus("Ready");
                return true;
            }
            UpdateStatus("Invalid Smagorinsky Constant");
            return false;
        }

        private static bool ValidateGrid(TextBox textBox)
        {
            try
            {
                if (string.IsNullOrWhiteSpace(textBox.Text) || !textBox.Text.Contains("x"))
                {
                    textBox.Text = $"{gridPointsX}x{gridPointsZ}";
                    UpdateStatus("Invalid Grid Format (e.g., 50x20)");
                    return false;
                }

                var parts = textBox.Text.Split('x');
                if (parts.Length != 2 || !int.TryParse(parts[0], out int x) || !int.TryParse(parts[1], out int z))
                {
                    textBox.Text = $"{gridPointsX}x{gridPointsZ}";
                    UpdateStatus("Invalid Grid Format (e.g., 50x20)");
                    return false;
                }

                int newX = Math.Max(20, Math.Min(100, x));
                int newZ = Math.Max(10, Math.Min(50, z));
                if (newX != gridPointsX || newZ != gridPointsZ)
                {
                    gridPointsX = newX;
                    gridPointsZ = newZ;
                    InitializeSimulation(); // Reinitialize arrays with new grid sizes
                }
                textBox.Text = $"{gridPointsX}x{gridPointsZ}";
                UpdateStatus("Ready");
                return true;
            }
            catch (Exception ex)
            {
                textBox.Text = $"{gridPointsX}x{gridPointsZ}";
                UpdateStatus($"Grid validation error: {ex.Message}");
                return false;
            }
        }

        private static bool ValidateInflow(TextBox textBox)
        {
            if (double.TryParse(textBox.Text, out double value))
            {
                riverInflow = Math.Max(0.01, Math.Min(100.0, value));
                textBox.Text = riverInflow.ToString("F2");
                UpdateStatus("Ready");
                return true;
            }
            UpdateStatus("Invalid River Inflow");
            return false;
        }

        private static bool ValidateTidal(TextBox textBox)
        {
            if (double.TryParse(textBox.Text, out double value))
            {
                tidalAmplitude = Math.Max(0.1, Math.Min(10.0, value));
                textBox.Text = tidalAmplitude.ToString("F2");
                UpdateStatus("Ready");
                return true;
            }
            UpdateStatus("Invalid Tidal Amplitude");
            return false;
        }

        private static bool ValidateAllInputs(TextBox smagorinskyTextBox, TextBox gridTextBox, TextBox inflowTextBox, TextBox tidalTextBox)
        {
            bool isValid = true;
            isValid &= ValidateSmagorinsky(smagorinskyTextBox);
            isValid &= ValidateGrid(gridTextBox);
            isValid &= ValidateInflow(inflowTextBox);
            isValid &= ValidateTidal(tidalTextBox);
            return isValid;
        }

        private static void UpdateStatus(string message)
        {
            if (statusLabel != null)
                statusLabel.Text = message;
        }

        private static void ResetSimulation()
        {
            simulationTime = 0.0;
            velocityU = new double[gridPointsX, gridPointsZ];
            velocityW = new double[gridPointsX, gridPointsZ];
            salinity = new double[gridPointsX, gridPointsZ];
            temperature = new double[gridPointsX, gridPointsZ];
            vorticity = new double[gridPointsX, gridPointsZ];
            eddyViscosity = new double[gridPointsX, gridPointsZ];

            // Initialize velocity (river inflow on left, tidal on right)
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    double x = i * dx;
                    velocityU[i, j] = (x < estuaryLength / 2) ? riverInflow / estuaryDepth : tidalAmplitude * Math.Cos(2 * Math.PI * simulationTime / tidalPeriod);
                    velocityW[i, j] = 0.0;
                    salinity[i, j] = (x / estuaryLength) * salinityOcean;
                    temperature[i, j] = 15.0 + (x / estuaryLength) * (temperatureOcean - 15.0);
                    vorticity[i, j] = 0.0;
                    eddyViscosity[i, j] = kinematicViscosity;
                    // Clamp initial values
                    velocityU[i, j] = Math.Max(-1.0, Math.Min(1.0, velocityU[i, j]));
                    velocityW[i, j] = Math.Max(-0.1, Math.Min(0.1, velocityW[i, j]));
                    salinity[i, j] = Math.Max(0.0, Math.Min(salinityOcean, salinity[i, j]));
                    temperature[i, j] = Math.Max(0.0, Math.Min(temperatureOcean, temperature[i, j]));
                }
            }
        }

        private static void UpdateSimulation()
        {
            if (salinity == null || temperature == null || salinity.GetLength(0) != gridPointsX || salinity.GetLength(1) != gridPointsZ)
            {
                UpdateStatus("Error: Arrays not initialized");
                MessageBox.Show("Simulation arrays not properly initialized", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return;
            }

            // Check CFL condition for stability
            double maxU = 0.0, maxW = 0.0, maxNu = kinematicViscosity;
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    maxU = Math.Max(maxU, Math.Abs(velocityU[i, j]));
                    maxW = Math.Max(maxW, Math.Abs(velocityW[i, j]));
                    maxNu = Math.Max(maxNu, eddyViscosity[i, j]);
                }
            }
            double cflAdv = Math.Min(dx / (maxU + 1e-10), dz / (maxW + 1e-10));
            double cflDiff = 0.5 * Math.Min(dx * dx, dz * dz) / (maxNu + 1e-10);
            double cfl = Math.Min(cflAdv, cflDiff);
            if (dt > cfl)
            {
                UpdateStatus($"Error: Time step too large (dt={dt:F2} > CFL={cfl:F2})");
                MessageBox.Show($"Time step too large: dt={dt:F2}s exceeds CFL limit {cfl:F2}s", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                simulationTimer.Stop();
                isSimulationRunning = false;
                return;
            }

            double[,] newVelocityU = new double[gridPointsX, gridPointsZ];
            double[,] newVelocityW = new double[gridPointsX, gridPointsZ];
            double[,] newSalinity = new double[gridPointsX, gridPointsZ];
            double[,] newTemperature = new double[gridPointsX, gridPointsZ];

            // Calculate density and baroclinic pressure gradient
            double[,] density = new double[gridPointsX, gridPointsZ];
            double[,] pressureGradientX = new double[gridPointsX, gridPointsZ];
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    density[i, j] = EqOfState.ComputeDensity(salinity[i, j], temperature[i, j], densityReference * g * (estuaryDepth - j * dz) / 1e4);
                    density[i, j] = Math.Max(1000.0, Math.Min(1030.0, density[i, j]));
                }
            }
            // Integrate pressure gradient over control volume (average density gradient)
            for (int i = 1; i < gridPointsX - 1; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    double rhoEast = 0.5 * (density[i, j] + density[i + 1, j]);
                    double rhoWest = 0.5 * (density[i - 1, j] + density[i, j]);
                    pressureGradientX[i, j] = g * (rhoEast - rhoWest) / dx;
                    pressureGradientX[i, j] = Math.Max(-100.0, Math.Min(100.0, pressureGradientX[i, j]));
                }
            }

            // Update eddy viscosity using Smagorinsky model
            for (int i = 1; i < gridPointsX - 1; i++)
            {
                for (int j = 1; j < gridPointsZ - 1; j++)
                {
                    double dudx = (velocityU[i + 1, j] - velocityU[i - 1, j]) / (2 * dx);
                    double dudz = (velocityU[i, j + 1] - velocityU[i, j - 1]) / (2 * dz);
                    double dwdx = (velocityW[i + 1, j] - velocityW[i - 1, j]) / (2 * dx);
                    double dwdz = (velocityW[i, j + 1] - velocityW[i, j - 1]) / (2 * dz);
                    double strainRate = Math.Sqrt(2 * (dudx * dudx + dwdz * dwdz) + (dudz + dwdx) * (dudz + dwdx));
                    double filterLength = Math.Sqrt(dx * dz);
                    eddyViscosity[i, j] = Math.Pow(smagorinskyConstant * filterLength, 2) * strainRate;
                    eddyViscosity[i, j] = Math.Max(kinematicViscosity, Math.Min(0.1, eddyViscosity[i, j]));
                }
            }

            // Finite Volume Method for momentum and scalar equations
            double cellArea = dx * dz; // Control volume area
            for (int i = 1; i < gridPointsX - 1; i++)
            {
                for (int j = 1; j < gridPointsZ - 1; j++)
                {
                    // Face velocities (interpolate to faces using upwind scheme)
                    double uEast = velocityU[i, j] >= 0 ? velocityU[i, j] : velocityU[i + 1, j];
                    double uWest = velocityU[i - 1, j] >= 0 ? velocityU[i - 1, j] : velocityU[i, j];
                    double wNorth = velocityW[i, j] >= 0 ? velocityW[i, j] : velocityW[i, j + 1];
                    double wSouth = velocityW[i, j - 1] >= 0 ? velocityW[i, j - 1] : velocityW[i, j];

                    // Face eddy viscosities (average)
                    double nuEast = 0.5 * (eddyViscosity[i, j] + eddyViscosity[i + 1, j]);
                    double nuWest = 0.5 * (eddyViscosity[i - 1, j] + eddyViscosity[i, j]);
                    double nuNorth = 0.5 * (eddyViscosity[i, j] + eddyViscosity[i, j + 1]);
                    double nuSouth = 0.5 * (eddyViscosity[i, j - 1] + eddyViscosity[i, j]);

                    // --- u-momentum equation ---
                    // Advective fluxes: ρu * u * A (A = dz for east/west faces)
                    double advFluxUEast = density[i, j] * uEast * velocityU[i, j] * dz;
                    double advFluxUWest = density[i - 1, j] * uWest * velocityU[i - 1, j] * dz;
                    double advFluxUNorth = density[i, j] * wNorth * velocityU[i, j] * dx;
                    double advFluxUSouth = density[i, j - 1] * wSouth * velocityU[i, j - 1] * dx;

                    // Diffusive fluxes: ν * (∂u/∂x) * A
                    double diffFluxUEast = nuEast * (velocityU[i + 1, j] - velocityU[i, j]) / dx * dz;
                    double diffFluxUWest = nuWest * (velocityU[i, j] - velocityU[i - 1, j]) / dx * dz;
                    double diffFluxUNorth = nuNorth * (velocityU[i, j + 1] - velocityU[i, j]) / dz * dx;
                    double diffFluxUSouth = nuSouth * (velocityU[i, j] - velocityU[i, j - 1]) / dz * dx;

                    // Net flux for u-momentum
                    double netFluxU = -(advFluxUEast - advFluxUWest + advFluxUNorth - advFluxUSouth) +
                                      (diffFluxUEast - diffFluxUWest + diffFluxUNorth - diffFluxUSouth);

                    // Source term: pressure gradient
                    double sourceU = -pressureGradientX[i, j] * cellArea / densityReference;

                    // Update u-velocity
                    newVelocityU[i, j] = velocityU[i, j] + dt * (netFluxU / (density[i, j] * cellArea) + sourceU / density[i, j]);

                    // --- w-momentum equation ---
                    // Advective fluxes: ρw * w * A
                    double advFluxWEast = density[i, j] * uEast * velocityW[i, j] * dz;
                    double advFluxWWest = density[i - 1, j] * uWest * velocityW[i - 1, j] * dz;
                    double advFluxWNorth = density[i, j] * wNorth * velocityW[i, j] * dx;
                    double advFluxWSouth = density[i, j - 1] * wSouth * velocityW[i, j - 1] * dx;

                    // Diffusive fluxes: ν * (∂w/∂x) * A
                    double diffFluxWEast = nuEast * (velocityW[i + 1, j] - velocityW[i, j]) / dx * dz;
                    double diffFluxWWest = nuWest * (velocityW[i, j] - velocityW[i - 1, j]) / dx * dz;
                    double diffFluxWNorth = nuNorth * (velocityW[i, j + 1] - velocityW[i, j]) / dz * dx;
                    double diffFluxWSouth = nuSouth * (velocityW[i, j] - velocityW[i, j - 1]) / dz * dx;

                    // Net flux for w-momentum
                    double netFluxW = -(advFluxWEast - advFluxWWest + advFluxWNorth - advFluxWSouth) +
                                      (diffFluxWEast - diffFluxWWest + diffFluxWNorth - diffFluxWSouth);

                    // Update w-velocity (no pressure gradient in z-direction)
                    newVelocityW[i, j] = velocityW[i, j] + dt * (netFluxW / (density[i, j] * cellArea));

                    // --- Salinity equation ---
                    // Advective fluxes: s * u * A
                    double sEast = velocityU[i, j] >= 0 ? salinity[i, j] : salinity[i + 1, j];
                    double sWest = velocityU[i - 1, j] >= 0 ? salinity[i - 1, j] : salinity[i, j];
                    double sNorth = velocityW[i, j] >= 0 ? salinity[i, j] : salinity[i, j + 1];
                    double sSouth = velocityW[i, j - 1] >= 0 ? salinity[i, j - 1] : salinity[i, j];

                    double advFluxSEast = sEast * uEast * dz;
                    double advFluxSWest = sWest * uWest * dz;
                    double advFluxSNorth = sNorth * wNorth * dx;
                    double advFluxSSouth = sSouth * wSouth * dx;

                    // Diffusive fluxes: ν * (∂s/∂x) * A
                    double diffFluxSEast = nuEast * (salinity[i + 1, j] - salinity[i, j]) / dx * dz;
                    double diffFluxSWest = nuWest * (salinity[i, j] - salinity[i - 1, j]) / dx * dz;
                    double diffFluxSNorth = nuNorth * (salinity[i, j + 1] - salinity[i, j]) / dz * dx;
                    double diffFluxSSouth = nuSouth * (salinity[i, j] - salinity[i, j - 1]) / dz * dx;

                    // Net flux for salinity
                    double netFluxS = -(advFluxSEast - advFluxSWest + advFluxSNorth - advFluxSSouth) +
                                      (diffFluxSEast - diffFluxSWest + diffFluxSNorth - diffFluxSSouth);

                    // Update salinity
                    newSalinity[i, j] = salinity[i, j] + dt * (netFluxS / cellArea);

                    // --- Temperature equation ---
                    // Advective fluxes: T * u * A
                    double tEast = velocityU[i, j] >= 0 ? temperature[i, j] : temperature[i + 1, j];
                    double tWest = velocityU[i - 1, j] >= 0 ? temperature[i - 1, j] : temperature[i, j];
                    double tNorth = velocityW[i, j] >= 0 ? temperature[i, j] : temperature[i, j + 1];
                    double tSouth = velocityW[i, j - 1] >= 0 ? temperature[i, j - 1] : temperature[i, j];

                    double advFluxTEast = tEast * uEast * dz;
                    double advFluxTWest = tWest * uWest * dz;
                    double advFluxTNorth = tNorth * wNorth * dx;
                    double advFluxTSouth = tSouth * wSouth * dx;

                    // Diffusive fluxes: ν * (∂T/∂x) * A
                    double diffFluxTEast = nuEast * (temperature[i + 1, j] - temperature[i, j]) / dx * dz;
                    double diffFluxTWest = nuWest * (temperature[i, j] - temperature[i - 1, j]) / dx * dz;
                    double diffFluxTNorth = nuNorth * (temperature[i, j + 1] - temperature[i, j]) / dz * dx;
                    double diffFluxTSouth = nuSouth * (temperature[i, j] - temperature[i, j - 1]) / dz * dx;

                    // Net flux for temperature
                    double netFluxT = -(advFluxTEast - advFluxTWest + advFluxTNorth - advFluxTSouth) +
                                      (diffFluxTEast - diffFluxTWest + diffFluxTNorth - diffFluxTSouth);

                    // Update temperature
                    newTemperature[i, j] = temperature[i, j] + dt * (netFluxT / cellArea);

                    // Clamp values
                    newVelocityU[i, j] = Math.Max(-1.0, Math.Min(1.0, newVelocityU[i, j]));
                    newVelocityW[i, j] = Math.Max(-0.1, Math.Min(0.1, newVelocityW[i, j]));
                    newSalinity[i, j] = Math.Max(0.0, Math.Min(salinityOcean, newSalinity[i, j]));
                    newTemperature[i, j] = Math.Max(0.0, Math.Min(temperatureOcean, newTemperature[i, j]));
                }
            }

            // Boundary conditions
            for (int j = 0; j < gridPointsZ; j++)
            {
                // West boundary (river inflow)
                newVelocityU[0, j] = riverInflow / estuaryDepth;
                newVelocityW[0, j] = 0.0;
                newSalinity[0, j] = 0.0;
                newTemperature[0, j] = 15.0;

                // East boundary (tidal forcing)
                newVelocityU[gridPointsX - 1, j] = tidalAmplitude * Math.Cos(2 * Math.PI * simulationTime / tidalPeriod);
                newVelocityW[gridPointsX - 1, j] = 0.0;
                newSalinity[gridPointsX - 1, j] = salinityOcean;
                newTemperature[gridPointsX - 1, j] = temperatureOcean;
            }
            for (int i = 0; i < gridPointsX; i++)
            {
                // Bottom boundary (no-slip)
                newVelocityW[i, 0] = 0.0;
                // Surface boundary (free-slip)
                newVelocityW[i, gridPointsZ - 1] = 0.0;
            }

            // Update vorticity: (∂w/∂x - ∂u/∂z)
            for (int i = 1; i < gridPointsX - 1; i++)
            {
                for (int j = 1; j < gridPointsZ - 1; j++)
                {
                    vorticity[i, j] = (velocityW[i + 1, j] - velocityW[i - 1, j]) / (2 * dx) -
                                      (velocityU[i, j + 1] - velocityU[i, j - 1]) / (2 * dz);
                    vorticity[i, j] = Math.Max(-10.0, Math.Min(10.0, vorticity[i, j]));
                }
            }

            // Copy new fields
            Array.Copy(newVelocityU, velocityU, gridPointsX * gridPointsZ);
            Array.Copy(newVelocityW, velocityW, gridPointsX * gridPointsZ);
            Array.Copy(newSalinity, salinity, gridPointsX * gridPointsZ);
            Array.Copy(newTemperature, temperature, gridPointsX * gridPointsZ);

            // Advance time
            simulationTime += dt;
        }

        private static void DrawSimulation(Graphics g, int width, int height)
        {
            g.Clear(Color.White);

            // Draw water background
            using (Brush waterBrush = new SolidBrush(Color.LightBlue))
            {
                g.FillRectangle(waterBrush, 0, 0, width, height);
            }

            // Draw vorticity (color map)
            double maxVorticity = 1.0; // Scale for visualization
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    int pixelX = i * width / gridPointsX;
                    int pixelY = height - (j + 1) * height / gridPointsZ; // Flip z-axis
                    pixelX = Math.Max(0, Math.Min(width - 1, pixelX));
                    pixelY = Math.Max(0, Math.Min(height - 1, pixelY));
                    int cellWidth = Math.Max(1, width / gridPointsX);
                    int cellHeight = Math.Max(1, height / gridPointsZ);
                    double vort = Math.Max(-maxVorticity, Math.Min(maxVorticity, vorticity[i, j]));
                    int colorValue = (int)(255 * (vort + maxVorticity) / (2 * maxVorticity));
                    colorValue = Math.Max(0, Math.Min(255, colorValue));
                    using (Brush vortBrush = new SolidBrush(Color.FromArgb(colorValue, 0, 255 - colorValue)))
                    {
                        g.FillRectangle(vortBrush, pixelX, pixelY, cellWidth, cellHeight);
                    }
                }
            }

            // Draw velocity vectors
            using (Pen velocityPen = new Pen(Color.Black, 1))
            {
                double maxVelocity = 0.5; // Scale for visualization
                for (int i = 0; i < gridPointsX; i += 2)
                {
                    for (int j = 0; j < gridPointsZ; j += 2)
                    {
                        int pixelX = i * width / gridPointsX + width / (2 * gridPointsX);
                        int pixelY = height - (j + 1) * height / gridPointsZ + height / (2 * gridPointsZ);
                        pixelX = Math.Max(0, Math.Min(width - 1, pixelX));
                        pixelY = Math.Max(0, Math.Min(height - 1, pixelY));
                        double u = Math.Max(-maxVelocity, Math.Min(maxVelocity, velocityU[i, j]));
                        double w = Math.Max(-maxVelocity, Math.Min(maxVelocity, velocityW[i, j]));
                        int arrowLengthX = (int)(u / maxVelocity * 10);
                        int arrowLengthY = (int)(w / maxVelocity * 10);
                        int endX = pixelX + arrowLengthX;
                        int endY = pixelY - arrowLengthY; // Negative due to flipped z-axis
                        endX = Math.Max(0, Math.Min(width - 1, endX));
                        endY = Math.Max(0, Math.Min(height - 1, endY));
                        g.DrawLine(velocityPen, pixelX, pixelY, endX, endY);
                        if (arrowLengthX != 0 || arrowLengthY != 0)
                        {
                            int arrowSize = 3;
                            double angle = Math.Atan2(-arrowLengthY, arrowLengthX);
                            int arrowX1 = endX - (int)(arrowSize * Math.Cos(angle + Math.PI / 6));
                            int arrowY1 = endY + (int)(arrowSize * Math.Sin(angle + Math.PI / 6));
                            int arrowX2 = endX - (int)(arrowSize * Math.Cos(angle - Math.PI / 6));
                            int arrowY2 = endY + (int)(arrowSize * Math.Sin(angle - Math.PI / 6));
                            arrowX1 = Math.Max(0, Math.Min(width - 1, arrowX1));
                            arrowY1 = Math.Max(0, Math.Min(height - 1, arrowY1));
                            arrowX2 = Math.Max(0, Math.Min(width - 1, arrowX2));
                            arrowY2 = Math.Max(0, Math.Min(height - 1, arrowY2));
                            g.DrawLine(velocityPen, endX, endY, arrowX1, arrowY1);
                            g.DrawLine(velocityPen, endX, endY, arrowX2, arrowY2);
                        }
                    }
                }
            }

            // Draw labels (using Segoe UI for mobile look)
            using (Font font = new Font("Segoe UI", 10F))
            {
                g.DrawString("River", font, Brushes.Black, 10, height - 20);
                g.DrawString("Ocean", font, Brushes.Black, width - 50, height - 20);
                g.DrawString("Vorticity (Blue-Red)", font, Brushes.Black, 10, 10);
                g.DrawString("Velocity Vectors (Black)", font, Brushes.Black, 10, 30);
                g.DrawString($"Time: {simulationTime:F1} s", font, Brushes.Black, width - 100, 10);
            }
        }
    }
}