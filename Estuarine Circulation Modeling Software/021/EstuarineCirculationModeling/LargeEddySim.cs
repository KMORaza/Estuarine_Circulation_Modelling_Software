using System;
using System.Drawing;
using System.Windows.Forms;

namespace EstuarineCirculationModeling
{
    public static class LargeEddySim
    {
        private static int gridPointsX = 50; // Grid points in x-direction
        private static int gridPointsZ = 20; // Grid points in z-direction
        private static readonly int maxTotalGridPoints = 4000; // Maximum total grid points (X * Z)
        private static double estuaryLength = 10000.0; // m
        private static double estuaryDepth = 10.0; // m
        private static double dx; // Spatial step x (m)
        private static double dz; // Spatial step z (m)
        private static double dt = 1.0; // Time step (s), will be adjusted dynamically
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
        private static string timeIntegrationScheme = "RK4"; // Time integration scheme (RK4 or CN)
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
            UpdateStatus($"Ready (dt={dt:F2}s)");

            // Show form
            lesForm.ShowDialog();
        }

        private static void InitializeSimulation()
        {
            // Ensure valid grid sizes
            if (gridPointsX <= 0 || gridPointsZ <= 0)
            {
                gridPointsX = 50; // Default fallback
                gridPointsZ = 20;
                UpdateStatus("Grid sizes reset to default (50x20)");
            }

            try
            {
                dx = estuaryLength / gridPointsX;
                dz = estuaryDepth / gridPointsZ;

                // Adjust dt based on CFL condition
                double maxU = riverInflow / estuaryDepth; // Initial estimate
                double maxW = 0.1; // Conservative estimate
                double maxNu = 0.1; // Max eddy viscosity from Smagorinsky
                double cflAdv = Math.Min(dx / (maxU + 1e-10), dz / (maxW + 1e-10));
                double cflDiff = 0.5 * Math.Min(dx * dx, dz * dz) / (maxNu + 1e-10);
                dt = 0.8 * (timeIntegrationScheme == "RK4" ? Math.Min(cflAdv, cflDiff) : cflAdv); // Safety factor 0.8

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
                AutoScroll = true,
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

            // Grid resolution X
            Label gridXLabel = new Label
            {
                Text = "Grid Points X (20-100):",
                Location = new Point(10, 70),
                AutoSize = true,
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black
            };
            TextBox gridXTextBox = new TextBox
            {
                Location = new Point(10, 95),
                Size = new Size(150, 25),
                Text = gridPointsX.ToString(),
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            gridXTextBox.TextChanged += (s, e) => ValidateGrid(gridXTextBox, null);

            // Grid resolution Z
            Label gridZLabel = new Label
            {
                Text = "Grid Points Z (10-50):",
                Location = new Point(170, 70),
                AutoSize = true,
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black
            };
            TextBox gridZTextBox = new TextBox
            {
                Location = new Point(170, 95),
                Size = new Size(150, 25),
                Text = gridPointsZ.ToString(),
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            gridZTextBox.TextChanged += (s, e) => ValidateGrid(null, gridZTextBox);

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

            // Time integration scheme
            Label schemeLabel = new Label
            {
                Text = "Time Integration:",
                Location = new Point(10, 250),
                AutoSize = true,
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black
            };
            ComboBox schemeComboBox = new ComboBox
            {
                Location = new Point(10, 275),
                Size = new Size(310, 25),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = new Font("Segoe UI", 10F),
                BackColor = Color.White,
                ForeColor = Color.Black
            };
            schemeComboBox.Items.AddRange(new string[] { "RK4", "CN" });
            schemeComboBox.SelectedIndex = 0; // Default to RK4
            schemeComboBox.SelectedIndexChanged += (s, e) =>
            {
                if (ValidateScheme(schemeComboBox))
                {
                    // Recalculate dt based on new scheme
                    double maxU = riverInflow / estuaryDepth;
                    double maxW = 0.1;
                    double maxNu = 0.1;
                    double cflAdv = Math.Min(dx / (maxU + 1e-10), dz / (maxW + 1e-10));
                    double cflDiff = 0.5 * Math.Min(dx * dx, dz * dz) / (maxNu + 1e-10);
                    dt = 0.8 * (timeIntegrationScheme == "RK4" ? Math.Min(cflAdv, cflDiff) : cflAdv);
                    UpdateStatus($"Time integration set to {timeIntegrationScheme} (dt={dt:F2}s)");
                }
            };

            // Apply button
            Button applyButton = new Button
            {
                Text = "Apply",
                Location = new Point(10, 310),
                Size = new Size(100, 30),
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = Color.FromArgb(200, 200, 200),
                ForeColor = Color.Black
            };
            applyButton.Click += (s, e) =>
            {
                if (ValidateAllInputs(smagorinskyTextBox, gridXTextBox, gridZTextBox, inflowTextBox, tidalTextBox, schemeComboBox))
                {
                    InitializeSimulation();
                    simulationPanel.Invalidate();
                    UpdateStatus($"Parameters applied (dt={dt:F2}s)");
                }
            };

            // Close button
            Button closeButton = new Button
            {
                Text = "Close",
                Location = new Point(120, 310),
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
                Location = new Point(10, 350),
                Size = new Size(100, 30),
                Font = new Font("Segoe UI", 12F, FontStyle.Bold),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = Color.FromArgb(200, 200, 200),
                ForeColor = Color.Black
            };
            startButton.Click += (s, e) =>
            {
                if (ValidateAllInputs(smagorinskyTextBox, gridXTextBox, gridZTextBox, inflowTextBox, tidalTextBox, schemeComboBox))
                {
                    if (salinity == null || temperature == null || salinity.GetLength(0) != gridPointsX || salinity.GetLength(1) != gridPointsZ)
                    {
                        InitializeSimulation();
                    }
                    simulationTimer.Start();
                    isSimulationRunning = true;
                    UpdateStatus($"Running (dt={dt:F2}s)");
                }
            };

            // Pause button
            Button pauseButton = new Button
            {
                Text = "Pause",
                Location = new Point(120, 350),
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
                UpdateStatus($"Paused (dt={dt:F2}s)");
            };

            // Reset button
            Button resetButton = new Button
            {
                Text = "Reset",
                Location = new Point(230, 350),
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
                schemeComboBox.SelectedIndex = 0; // Reset to RK4
                timeIntegrationScheme = "RK4";
                gridXTextBox.Text = gridPointsX.ToString();
                gridZTextBox.Text = gridPointsZ.ToString();
                simulationPanel.Invalidate();
                UpdateStatus($"Reset (dt={dt:F2}s)");
            };

            // Status label
            statusLabel = new Label
            {
                Text = $"Ready (dt={dt:F2}s)",
                Location = new Point(10, 390),
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
                UpdateStatus($"Running (Time: {simulationTime:F1} s, Scheme: {timeIntegrationScheme}, dt={dt:F2}s)");
            };

            // Add controls to control panel
            controlPanel.Controls.AddRange(new Control[] { smagorinskyLabel, smagorinskyTextBox, gridXLabel, gridXTextBox, gridZLabel, gridZTextBox, inflowLabel, inflowTextBox, tidalLabel, tidalTextBox, schemeLabel, schemeComboBox, applyButton, closeButton, startButton, pauseButton, resetButton, statusLabel });

            // Add panels to form
            lesForm.Controls.Add(borderPanel);
        }

        private static bool ValidateSmagorinsky(TextBox textBox)
        {
            if (double.TryParse(textBox.Text, out double value))
            {
                smagorinskyConstant = Math.Max(0.01, Math.Min(0.5, value));
                textBox.Text = smagorinskyConstant.ToString("F2");
                UpdateStatus($"Ready (dt={dt:F2}s)");
                return true;
            }
            UpdateStatus("Invalid Smagorinsky Constant");
            textBox.Text = smagorinskyConstant.ToString("F2");
            return false;
        }

        private static bool ValidateGrid(TextBox gridXTextBox, TextBox gridZTextBox)
        {
            try
            {
                bool isValid = true;
                int newX = gridPointsX;
                int newZ = gridPointsZ;

                if (gridXTextBox != null)
                {
                    if (int.TryParse(gridXTextBox.Text, out int x))
                    {
                        newX = Math.Max(20, Math.Min(100, x));
                        gridXTextBox.Text = newX.ToString();
                    }
                    else
                    {
                        isValid = false;
                        gridXTextBox.Text = gridPointsX.ToString();
                        UpdateStatus("Invalid Grid Points X (20-100)");
                    }
                }

                if (gridZTextBox != null)
                {
                    if (int.TryParse(gridZTextBox.Text, out int z))
                    {
                        newZ = Math.Max(10, Math.Min(50, z));
                        gridZTextBox.Text = newZ.ToString();
                    }
                    else
                    {
                        isValid = false;
                        gridZTextBox.Text = gridPointsZ.ToString();
                        UpdateStatus("Invalid Grid Points Z (10-50)");
                    }
                }

                // Check total grid points limit
                if (newX * newZ > maxTotalGridPoints)
                {
                    isValid = false;
                    // Adjust gridZ first to minimize total points while respecting minimums
                    newZ = Math.Max(10, Math.Min(50, maxTotalGridPoints / newX));
                    if (newX * newZ > maxTotalGridPoints)
                    {
                        // If still over, adjust gridX
                        newX = Math.Max(20, Math.Min(100, maxTotalGridPoints / newZ));
                        newZ = Math.Max(10, Math.Min(50, maxTotalGridPoints / newX));
                    }
                    if (gridXTextBox != null) gridXTextBox.Text = newX.ToString();
                    if (gridZTextBox != null) gridZTextBox.Text = newZ.ToString();
                    UpdateStatus($"Grid adjusted to {newX}x{newZ} (total points ≤ {maxTotalGridPoints})");
                }

                if (newX != gridPointsX || newZ != gridPointsZ)
                {
                    gridPointsX = newX;
                    gridPointsZ = newZ;
                    InitializeSimulation(); // Reinitialize arrays and adjust dt
                    UpdateStatus($"Grid updated to {gridPointsX}x{gridPointsZ} (dt={dt:F2}s)");
                }

                return isValid;
            }
            catch (Exception ex)
            {
                if (gridXTextBox != null) gridXTextBox.Text = gridPointsX.ToString();
                if (gridZTextBox != null) gridZTextBox.Text = gridPointsZ.ToString();
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
                UpdateStatus($"Ready (dt={dt:F2}s)");
                return true;
            }
            UpdateStatus("Invalid River Inflow");
            textBox.Text = riverInflow.ToString("F2");
            return false;
        }

        private static bool ValidateTidal(TextBox textBox)
        {
            if (double.TryParse(textBox.Text, out double value))
            {
                tidalAmplitude = Math.Max(0.1, Math.Min(10.0, value));
                textBox.Text = tidalAmplitude.ToString("F2");
                UpdateStatus($"Ready (dt={dt:F2}s)");
                return true;
            }
            UpdateStatus("Invalid Tidal Amplitude");
            textBox.Text = tidalAmplitude.ToString("F2");
            return false;
        }

        private static bool ValidateScheme(ComboBox comboBox)
        {
            if (comboBox.SelectedItem == null)
            {
                comboBox.SelectedIndex = 0;
                timeIntegrationScheme = "RK4";
                UpdateStatus($"Time integration scheme reset to RK4 (dt={dt:F2}s)");
                return false;
            }
            timeIntegrationScheme = comboBox.SelectedItem.ToString();
            UpdateStatus($"Ready (dt={dt:F2}s)");
            return true;
        }

        private static bool ValidateAllInputs(TextBox smagorinskyTextBox, TextBox gridXTextBox, TextBox gridZTextBox, TextBox inflowTextBox, TextBox tidalTextBox, ComboBox schemeComboBox)
        {
            bool isValid = true;
            isValid &= ValidateSmagorinsky(smagorinskyTextBox);
            isValid &= ValidateGrid(gridXTextBox, gridZTextBox);
            isValid &= ValidateInflow(inflowTextBox);
            isValid &= ValidateTidal(tidalTextBox);
            isValid &= ValidateScheme(schemeComboBox);
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
            timeIntegrationScheme = "RK4"; // Reset to RK4
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
                    velocityU[i, j] = Math.Max(-1.0, Math.Min(1.0, velocityU[i, j]));
                    velocityW[i, j] = Math.Max(-0.1, Math.Min(0.1, velocityW[i, j]));
                    salinity[i, j] = Math.Max(0.0, Math.Min(salinityOcean, salinity[i, j]));
                    temperature[i, j] = Math.Max(0.0, Math.Min(temperatureOcean, temperature[i, j]));
                }
            }
        }

        private static void ComputeTendencies(double[,] u, double[,] w, double[,] s, double[,] t, double[,] ev, double[,] density, double[,] pressureGradientX,
                                             out double[,] du_dt, out double[,] dw_dt, out double[,] ds_dt, out double[,] dt_dt)
        {
            du_dt = new double[gridPointsX, gridPointsZ];
            dw_dt = new double[gridPointsX, gridPointsZ];
            ds_dt = new double[gridPointsX, gridPointsZ];
            dt_dt = new double[gridPointsX, gridPointsZ];
            double cellArea = dx * dz;

            for (int i = 1; i < gridPointsX - 1; i++)
            {
                for (int j = 1; j < gridPointsZ - 1; j++)
                {
                    double uEast = u[i, j] >= 0 ? u[i, j] : u[i + 1, j];
                    double uWest = u[i - 1, j] >= 0 ? u[i - 1, j] : u[i, j];
                    double wNorth = w[i, j] >= 0 ? w[i, j] : w[i, j + 1];
                    double wSouth = w[i, j - 1] >= 0 ? w[i, j - 1] : w[i, j];

                    double nuEast = 0.5 * (ev[i, j] + ev[i + 1, j]);
                    double nuWest = 0.5 * (ev[i - 1, j] + ev[i, j]);
                    double nuNorth = 0.5 * (ev[i, j] + ev[i, j + 1]);
                    double nuSouth = 0.5 * (ev[i, j - 1] + ev[i, j]);

                    double advFluxUEast = density[i, j] * uEast * u[i, j] * dz;
                    double advFluxUWest = density[i - 1, j] * uWest * u[i - 1, j] * dz;
                    double advFluxUNorth = density[i, j] * wNorth * u[i, j] * dx;
                    double advFluxUSouth = density[i, j - 1] * wSouth * u[i, j - 1] * dx;
                    double diffFluxUEast = nuEast * (u[i + 1, j] - u[i, j]) / dx * dz;
                    double diffFluxUWest = nuWest * (u[i, j] - u[i - 1, j]) / dx * dz;
                    double diffFluxUNorth = nuNorth * (u[i, j + 1] - u[i, j]) / dz * dx;
                    double diffFluxUSouth = nuSouth * (u[i, j] - u[i, j - 1]) / dz * dx;
                    double netFluxU = -(advFluxUEast - advFluxUWest + advFluxUNorth - advFluxUSouth) +
                                      (diffFluxUEast - diffFluxUWest + diffFluxUNorth - diffFluxUSouth);
                    double sourceU = -pressureGradientX[i, j] * cellArea / densityReference;
                    du_dt[i, j] = (netFluxU / (density[i, j] * cellArea) + sourceU / density[i, j]);

                    double advFluxWEast = density[i, j] * uEast * w[i, j] * dz;
                    double advFluxWWest = density[i - 1, j] * uWest * w[i - 1, j] * dz;
                    double advFluxWNorth = density[i, j] * wNorth * w[i, j] * dx;
                    double advFluxWSouth = density[i, j - 1] * wSouth * w[i, j - 1] * dx;
                    double diffFluxWEast = nuEast * (w[i + 1, j] - w[i, j]) / dx * dz;
                    double diffFluxWWest = nuWest * (w[i, j] - w[i - 1, j]) / dx * dz;
                    double diffFluxWNorth = nuNorth * (w[i, j + 1] - w[i, j]) / dz * dx;
                    double diffFluxWSouth = nuSouth * (w[i, j] - w[i, j - 1]) / dz * dx;
                    double netFluxW = -(advFluxWEast - advFluxWWest + advFluxWNorth - advFluxWSouth) +
                                      (diffFluxWEast - diffFluxWWest + diffFluxWNorth - diffFluxWSouth);
                    dw_dt[i, j] = (netFluxW / (density[i, j] * cellArea));

                    double sEast = u[i, j] >= 0 ? s[i, j] : s[i + 1, j];
                    double sWest = u[i - 1, j] >= 0 ? s[i - 1, j] : s[i, j];
                    double sNorth = w[i, j] >= 0 ? s[i, j] : s[i, j + 1];
                    double sSouth = w[i, j - 1] >= 0 ? s[i, j - 1] : s[i, j];
                    double advFluxSEast = sEast * uEast * dz;
                    double advFluxSWest = sWest * uWest * dz;
                    double advFluxSNorth = sNorth * wNorth * dx;
                    double advFluxSSouth = sSouth * wSouth * dx;
                    double diffFluxSEast = nuEast * (s[i + 1, j] - s[i, j]) / dx * dz;
                    double diffFluxSWest = nuWest * (s[i, j] - s[i - 1, j]) / dx * dz;
                    double diffFluxSNorth = nuNorth * (s[i, j + 1] - s[i, j]) / dz * dx;
                    double diffFluxSSouth = nuSouth * (s[i, j] - s[i, j - 1]) / dz * dx;
                    double netFluxS = -(advFluxSEast - advFluxSWest + advFluxSNorth - advFluxSSouth) +
                                      (diffFluxSEast - diffFluxSWest + diffFluxSNorth - diffFluxSSouth);
                    ds_dt[i, j] = (netFluxS / cellArea);

                    double tEast = u[i, j] >= 0 ? t[i, j] : t[i + 1, j];
                    double tWest = u[i - 1, j] >= 0 ? t[i - 1, j] : t[i, j];
                    double tNorth = w[i, j] >= 0 ? t[i, j] : t[i, j + 1];
                    double tSouth = w[i, j - 1] >= 0 ? t[i, j - 1] : t[i, j];
                    double advFluxTEast = tEast * uEast * dz;
                    double advFluxTWest = tWest * uWest * dz;
                    double advFluxTNorth = tNorth * wNorth * dx;
                    double advFluxTSouth = tSouth * wSouth * dx;
                    double diffFluxTEast = nuEast * (t[i + 1, j] - t[i, j]) / dx * dz;
                    double diffFluxTWest = nuWest * (t[i, j] - t[i - 1, j]) / dx * dz;
                    double diffFluxTNorth = nuNorth * (t[i, j + 1] - t[i, j]) / dz * dx;
                    double diffFluxTSouth = nuSouth * (t[i, j] - t[i, j - 1]) / dz * dx;
                    double netFluxT = -(advFluxTEast - advFluxTWest + advFluxTNorth - advFluxTSouth) +
                                      (diffFluxTEast - diffFluxTWest + diffFluxTNorth - diffFluxTSouth);
                    dt_dt[i, j] = (netFluxT / cellArea);
                }
            }
        }

        private static void SolveTridiagonal(double[] a, double[] b, double[] c, double[] d, double[] x, int n)
        {
            double[] cPrime = new double[n];
            double[] dPrime = new double[n];
            cPrime[0] = c[0] / b[0];
            dPrime[0] = d[0] / b[0];
            for (int i = 1; i < n; i++)
            {
                double m = b[i] - a[i] * cPrime[i - 1];
                cPrime[i] = c[i] / m;
                dPrime[i] = (d[i] - a[i] * dPrime[i - 1]) / m;
            }
            x[n - 1] = dPrime[n - 1];
            for (int i = n - 2; i >= 0; i--)
            {
                x[i] = dPrime[i] - cPrime[i] * x[i + 1];
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

            // Check CFL condition and adjust dt
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
            double cfl = timeIntegrationScheme == "RK4" ? Math.Min(cflAdv, cflDiff) : cflAdv;
            if (dt > cfl)
            {
                dt = 0.8 * cfl; // Adjust dt with safety factor
                UpdateStatus($"Time step adjusted to dt={dt:F2}s for stability");
            }

            // Calculate density and pressure gradient
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

            double[,] newVelocityU = new double[gridPointsX, gridPointsZ];
            double[,] newVelocityW = new double[gridPointsX, gridPointsZ];
            double[,] newSalinity = new double[gridPointsX, gridPointsZ];
            double[,] newTemperature = new double[gridPointsX, gridPointsZ];

            if (timeIntegrationScheme == "RK4")
            {
                double[,] u_k1 = velocityU, w_k1 = velocityW, s_k1 = salinity, t_k1 = temperature;
                double[,] du_dt1, dw_dt1, ds_dt1, dt_dt1;
                ComputeTendencies(u_k1, w_k1, s_k1, t_k1, eddyViscosity, density, pressureGradientX, out du_dt1, out dw_dt1, out ds_dt1, out dt_dt1);

                double[,] u_k2 = new double[gridPointsX, gridPointsZ];
                double[,] w_k2 = new double[gridPointsX, gridPointsZ];
                double[,] s_k2 = new double[gridPointsX, gridPointsZ];
                double[,] t_k2 = new double[gridPointsX, gridPointsZ];
                for (int i = 0; i < gridPointsX; i++)
                {
                    for (int j = 0; j < gridPointsZ; j++)
                    {
                        u_k2[i, j] = velocityU[i, j] + 0.5 * dt * du_dt1[i, j];
                        w_k2[i, j] = velocityW[i, j] + 0.5 * dt * dw_dt1[i, j];
                        s_k2[i, j] = salinity[i, j] + 0.5 * dt * ds_dt1[i, j];
                        t_k2[i, j] = temperature[i, j] + 0.5 * dt * dt_dt1[i, j];
                    }
                }
                double[,] du_dt2, dw_dt2, ds_dt2, dt_dt2;
                ComputeTendencies(u_k2, w_k2, s_k2, t_k2, eddyViscosity, density, pressureGradientX, out du_dt2, out dw_dt2, out ds_dt2, out dt_dt2);

                double[,] u_k3 = new double[gridPointsX, gridPointsZ];
                double[,] w_k3 = new double[gridPointsX, gridPointsZ];
                double[,] s_k3 = new double[gridPointsX, gridPointsZ];
                double[,] t_k3 = new double[gridPointsX, gridPointsZ];
                for (int i = 0; i < gridPointsX; i++)
                {
                    for (int j = 0; j < gridPointsZ; j++)
                    {
                        u_k3[i, j] = velocityU[i, j] + 0.5 * dt * du_dt2[i, j];
                        w_k3[i, j] = velocityW[i, j] + 0.5 * dt * dw_dt2[i, j];
                        s_k3[i, j] = salinity[i, j] + 0.5 * dt * ds_dt2[i, j];
                        t_k3[i, j] = temperature[i, j] + 0.5 * dt * dt_dt2[i, j];
                    }
                }
                double[,] du_dt3, dw_dt3, ds_dt3, dt_dt3;
                ComputeTendencies(u_k3, w_k3, s_k3, t_k3, eddyViscosity, density, pressureGradientX, out du_dt3, out dw_dt3, out ds_dt3, out dt_dt3);

                double[,] u_k4 = new double[gridPointsX, gridPointsZ];
                double[,] w_k4 = new double[gridPointsX, gridPointsZ];
                double[,] s_k4 = new double[gridPointsX, gridPointsZ];
                double[,] t_k4 = new double[gridPointsX, gridPointsZ];
                for (int i = 0; i < gridPointsX; i++)
                {
                    for (int j = 0; j < gridPointsZ; j++)
                    {
                        u_k4[i, j] = velocityU[i, j] + dt * du_dt3[i, j];
                        w_k4[i, j] = velocityW[i, j] + dt * dw_dt3[i, j];
                        s_k4[i, j] = salinity[i, j] + dt * ds_dt3[i, j];
                        t_k4[i, j] = temperature[i, j] + dt * dt_dt3[i, j];
                    }
                }
                double[,] du_dt4, dw_dt4, ds_dt4, dt_dt4;
                ComputeTendencies(u_k4, w_k4, s_k4, t_k4, eddyViscosity, density, pressureGradientX, out du_dt4, out dw_dt4, out ds_dt4, out dt_dt4);

                for (int i = 1; i < gridPointsX - 1; i++)
                {
                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        newVelocityU[i, j] = velocityU[i, j] + (dt / 6.0) * (du_dt1[i, j] + 2 * du_dt2[i, j] + 2 * du_dt3[i, j] + du_dt4[i, j]);
                        newVelocityW[i, j] = velocityW[i, j] + (dt / 6.0) * (dw_dt1[i, j] + 2 * dw_dt2[i, j] + 2 * dw_dt3[i, j] + dw_dt4[i, j]);
                        newSalinity[i, j] = salinity[i, j] + (dt / 6.0) * (ds_dt1[i, j] + 2 * ds_dt2[i, j] + 2 * ds_dt3[i, j] + ds_dt4[i, j]);
                        newTemperature[i, j] = temperature[i, j] + (dt / 6.0) * (dt_dt1[i, j] + 2 * dt_dt2[i, j] + 2 * dt_dt3[i, j] + dt_dt4[i, j]);
                    }
                }
            }
            else // CN
            {
                double[,] tempU = new double[gridPointsX, gridPointsZ];
                double[,] tempW = new double[gridPointsX, gridPointsZ];
                double[,] tempS = new double[gridPointsX, gridPointsZ];
                double[,] tempT = new double[gridPointsX, gridPointsZ];
                double[,] du_dt, dw_dt, ds_dt, dt_dt;
                ComputeTendencies(velocityU, velocityW, salinity, temperature, eddyViscosity, density, pressureGradientX, out du_dt, out dw_dt, out ds_dt, out dt_dt);

                for (int i = 1; i < gridPointsX - 1; i++)
                {
                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        double cellArea = dx * dz;
                        double advFluxUEast = density[i, j] * (velocityU[i, j] >= 0 ? velocityU[i, j] : velocityU[i + 1, j]) * velocityU[i, j] * dz;
                        double advFluxUWest = density[i - 1, j] * (velocityU[i - 1, j] >= 0 ? velocityU[i - 1, j] : velocityU[i, j]) * velocityU[i - 1, j] * dz;
                        double advFluxUNorth = density[i, j] * (velocityW[i, j] >= 0 ? velocityW[i, j] : velocityW[i, j + 1]) * velocityU[i, j] * dx;
                        double advFluxUSouth = density[i, j - 1] * (velocityW[i, j - 1] >= 0 ? velocityW[i, j - 1] : velocityW[i, j]) * velocityU[i, j - 1] * dx;
                        double netAdvFluxU = -(advFluxUEast - advFluxUWest + advFluxUNorth - advFluxUSouth);
                        double sourceU = -pressureGradientX[i, j] * cellArea / densityReference;
                        tempU[i, j] = velocityU[i, j] + dt * (netAdvFluxU / (density[i, j] * cellArea) + sourceU / density[i, j]);

                        double advFluxWEast = density[i, j] * (velocityU[i, j] >= 0 ? velocityU[i, j] : velocityU[i + 1, j]) * velocityW[i, j] * dz;
                        double advFluxWWest = density[i - 1, j] * (velocityU[i - 1, j] >= 0 ? velocityU[i - 1, j] : velocityU[i, j]) * velocityW[i - 1, j] * dz;
                        double advFluxWNorth = density[i, j] * (velocityW[i, j] >= 0 ? velocityW[i, j] : velocityW[i, j + 1]) * velocityW[i, j] * dx;
                        double advFluxWSouth = density[i, j - 1] * (velocityW[i, j - 1] >= 0 ? velocityW[i, j - 1] : velocityW[i, j]) * velocityW[i, j - 1] * dx;
                        double netAdvFluxW = -(advFluxWEast - advFluxWWest + advFluxWNorth - advFluxWSouth);
                        tempW[i, j] = velocityW[i, j] + dt * (netAdvFluxW / (density[i, j] * cellArea));

                        double sEast = velocityU[i, j] >= 0 ? salinity[i, j] : salinity[i + 1, j];
                        double sWest = velocityU[i - 1, j] >= 0 ? salinity[i - 1, j] : salinity[i, j];
                        double sNorth = velocityW[i, j] >= 0 ? salinity[i, j] : salinity[i, j + 1];
                        double sSouth = velocityW[i, j - 1] >= 0 ? salinity[i, j - 1] : salinity[i, j];
                        double advFluxSEast = sEast * (velocityU[i, j] >= 0 ? velocityU[i, j] : velocityU[i + 1, j]) * dz;
                        double advFluxSWest = sWest * (velocityU[i - 1, j] >= 0 ? velocityU[i - 1, j] : velocityU[i, j]) * dz;
                        double advFluxSNorth = sNorth * (velocityW[i, j] >= 0 ? velocityW[i, j] : velocityW[i, j + 1]) * dx;
                        double advFluxSSouth = sSouth * (velocityW[i, j - 1] >= 0 ? velocityW[i, j - 1] : velocityW[i, j]) * dx;
                        double netAdvFluxS = -(advFluxSEast - advFluxSWest + advFluxSNorth - advFluxSSouth);
                        tempS[i, j] = salinity[i, j] + dt * (netAdvFluxS / cellArea);

                        double tEast = velocityU[i, j] >= 0 ? temperature[i, j] : temperature[i + 1, j];
                        double tWest = velocityU[i - 1, j] >= 0 ? temperature[i - 1, j] : temperature[i, j];
                        double tNorth = velocityW[i, j] >= 0 ? temperature[i, j] : temperature[i, j + 1];
                        double tSouth = velocityW[i, j - 1] >= 0 ? temperature[i, j - 1] : temperature[i, j];
                        double advFluxTEast = tEast * (velocityU[i, j] >= 0 ? velocityU[i, j] : velocityU[i + 1, j]) * dz;
                        double advFluxTWest = tWest * (velocityU[i - 1, j] >= 0 ? velocityU[i - 1, j] : velocityU[i, j]) * dz;
                        double advFluxTNorth = tNorth * (velocityW[i, j] >= 0 ? velocityW[i, j] : velocityW[i, j + 1]) * dx;
                        double advFluxTSouth = tSouth * (velocityW[i, j - 1] >= 0 ? velocityW[i, j - 1] : velocityW[i, j]) * dx;
                        double netAdvFluxT = -(advFluxTEast - advFluxTWest + advFluxTNorth - advFluxTSouth);
                        tempT[i, j] = temperature[i, j] + dt * (netAdvFluxT / cellArea);
                    }
                }

                for (int j = 1; j < gridPointsZ - 1; j++)
                {
                    double[] a = new double[gridPointsX - 2];
                    double[] b = new double[gridPointsX - 2];
                    double[] c = new double[gridPointsX - 2];
                    double[] d = new double[gridPointsX - 2];
                    double[] x = new double[gridPointsX - 2];

                    for (int i = 1; i < gridPointsX - 1; i++)
                    {
                        double nuEast = 0.5 * (eddyViscosity[i, j] + eddyViscosity[i + 1, j]);
                        double nuWest = 0.5 * (eddyViscosity[i - 1, j] + eddyViscosity[i, j]);
                        double alpha = 0.5 * dt * dz / (dx * dx);
                        a[i - 1] = -alpha * nuWest;
                        b[i - 1] = 1.0 + alpha * (nuWest + nuEast);
                        c[i - 1] = -alpha * nuEast;
                        d[i - 1] = tempU[i, j];
                    }
                    SolveTridiagonal(a, b, c, d, x, gridPointsX - 2);
                    for (int i = 1; i < gridPointsX - 1; i++)
                        newVelocityU[i, j] = x[i - 1];

                    for (int i = 1; i < gridPointsX - 1; i++)
                    {
                        double nuEast = 0.5 * (eddyViscosity[i, j] + eddyViscosity[i + 1, j]);
                        double nuWest = 0.5 * (eddyViscosity[i - 1, j] + eddyViscosity[i, j]);
                        double alpha = 0.5 * dt * dz / (dx * dx);
                        a[i - 1] = -alpha * nuWest;
                        b[i - 1] = 1.0 + alpha * (nuWest + nuEast);
                        c[i - 1] = -alpha * nuEast;
                        d[i - 1] = tempW[i, j];
                    }
                    SolveTridiagonal(a, b, c, d, x, gridPointsX - 2);
                    for (int i = 1; i < gridPointsX - 1; i++)
                        newVelocityW[i, j] = x[i - 1];

                    for (int i = 1; i < gridPointsX - 1; i++)
                    {
                        double nuEast = 0.5 * (eddyViscosity[i, j] + eddyViscosity[i + 1, j]);
                        double nuWest = 0.5 * (eddyViscosity[i - 1, j] + eddyViscosity[i, j]);
                        double alpha = 0.5 * dt * dz / (dx * dx);
                        a[i - 1] = -alpha * nuWest;
                        b[i - 1] = 1.0 + alpha * (nuWest + nuEast);
                        c[i - 1] = -alpha * nuEast;
                        d[i - 1] = tempS[i, j];
                    }
                    SolveTridiagonal(a, b, c, d, x, gridPointsX - 2);
                    for (int i = 1; i < gridPointsX - 1; i++)
                        newSalinity[i, j] = x[i - 1];

                    for (int i = 1; i < gridPointsX - 1; i++)
                    {
                        double nuEast = 0.5 * (eddyViscosity[i, j] + eddyViscosity[i + 1, j]);
                        double nuWest = 0.5 * (eddyViscosity[i - 1, j] + eddyViscosity[i, j]);
                        double alpha = 0.5 * dt * dz / (dx * dx);
                        a[i - 1] = -alpha * nuWest;
                        b[i - 1] = 1.0 + alpha * (nuWest + nuEast);
                        c[i - 1] = -alpha * nuEast;
                        d[i - 1] = tempT[i, j];
                    }
                    SolveTridiagonal(a, b, c, d, x, gridPointsX - 2);
                    for (int i = 1; i < gridPointsX - 1; i++)
                        newTemperature[i, j] = x[i - 1];
                }

                for (int i = 1; i < gridPointsX - 1; i++)
                {
                    double[] a = new double[gridPointsZ - 2];
                    double[] b = new double[gridPointsZ - 2];
                    double[] c = new double[gridPointsZ - 2];
                    double[] d = new double[gridPointsZ - 2];
                    double[] x = new double[gridPointsZ - 2];

                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        double nuNorth = 0.5 * (eddyViscosity[i, j] + eddyViscosity[i, j + 1]);
                        double nuSouth = 0.5 * (eddyViscosity[i, j - 1] + eddyViscosity[i, j]);
                        double alpha = 0.5 * dt * dx / (dz * dz);
                        a[j - 1] = -alpha * nuSouth;
                        b[j - 1] = 1.0 + alpha * (nuSouth + nuNorth);
                        c[j - 1] = -alpha * nuNorth;
                        d[j - 1] = newVelocityU[i, j];
                    }
                    SolveTridiagonal(a, b, c, d, x, gridPointsZ - 2);
                    for (int j = 1; j < gridPointsZ - 1; j++)
                        newVelocityU[i, j] = x[j - 1];

                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        double nuNorth = 0.5 * (eddyViscosity[i, j] + eddyViscosity[i, j + 1]);
                        double nuSouth = 0.5 * (eddyViscosity[i, j - 1] + eddyViscosity[i, j]);
                        double alpha = 0.5 * dt * dx / (dz * dz);
                        a[j - 1] = -alpha * nuSouth;
                        b[j - 1] = 1.0 + alpha * (nuSouth + nuNorth);
                        c[j - 1] = -alpha * nuNorth;
                        d[j - 1] = newVelocityW[i, j];
                    }
                    SolveTridiagonal(a, b, c, d, x, gridPointsZ - 2);
                    for (int j = 1; j < gridPointsZ - 1; j++)
                        newVelocityW[i, j] = x[j - 1];

                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        double nuNorth = 0.5 * (eddyViscosity[i, j] + eddyViscosity[i, j + 1]);
                        double nuSouth = 0.5 * (eddyViscosity[i, j - 1] + eddyViscosity[i, j]);
                        double alpha = 0.5 * dt * dx / (dz * dz);
                        a[j - 1] = -alpha * nuSouth;
                        b[j - 1] = 1.0 + alpha * (nuSouth + nuNorth);
                        c[j - 1] = -alpha * nuNorth;
                        d[j - 1] = newSalinity[i, j];
                    }
                    SolveTridiagonal(a, b, c, d, x, gridPointsZ - 2);
                    for (int j = 1; j < gridPointsZ - 1; j++)
                        newSalinity[i, j] = x[j - 1];

                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        double nuNorth = 0.5 * (eddyViscosity[i, j] + eddyViscosity[i, j + 1]);
                        double nuSouth = 0.5 * (eddyViscosity[i, j - 1] + eddyViscosity[i, j]);
                        double alpha = 0.5 * dt * dx / (dz * dz);
                        a[j - 1] = -alpha * nuSouth;
                        b[j - 1] = 1.0 + alpha * (nuSouth + nuNorth);
                        c[j - 1] = -alpha * nuNorth;
                        d[j - 1] = newTemperature[i, j];
                    }
                    SolveTridiagonal(a, b, c, d, x, gridPointsZ - 2);
                    for (int j = 1; j < gridPointsZ - 1; j++)
                        newTemperature[i, j] = x[j - 1];
                }
            }

            // Boundary conditions
            for (int j = 0; j < gridPointsZ; j++)
            {
                newVelocityU[0, j] = riverInflow / estuaryDepth;
                newVelocityW[0, j] = 0.0;
                newSalinity[0, j] = 0.0;
                newTemperature[0, j] = 15.0;
                newVelocityU[gridPointsX - 1, j] = tidalAmplitude * Math.Cos(2 * Math.PI * simulationTime / tidalPeriod);
                newVelocityW[gridPointsX - 1, j] = 0.0;
                newSalinity[gridPointsX - 1, j] = salinityOcean;
                newTemperature[gridPointsX - 1, j] = temperatureOcean;
            }
            for (int i = 0; i < gridPointsX; i++)
            {
                newVelocityW[i, 0] = 0.0; // Bottom
                newVelocityW[i, gridPointsZ - 1] = 0.0; // Surface
            }

            // Update vorticity
            for (int i = 1; i < gridPointsX - 1; i++)
            {
                for (int j = 1; j < gridPointsZ - 1; j++)
                {
                    vorticity[i, j] = (newVelocityW[i + 1, j] - newVelocityW[i - 1, j]) / (2 * dx) -
                                      (newVelocityU[i, j + 1] - newVelocityU[i, j - 1]) / (2 * dz);
                    vorticity[i, j] = Math.Max(-10.0, Math.Min(10.0, vorticity[i, j]));
                }
            }

            // Clamp values
            for (int i = 1; i < gridPointsX - 1; i++)
            {
                for (int j = 1; j < gridPointsZ - 1; j++)
                {
                    newVelocityU[i, j] = Math.Max(-1.0, Math.Min(1.0, newVelocityU[i, j]));
                    newVelocityW[i, j] = Math.Max(-0.1, Math.Min(0.1, newVelocityW[i, j]));
                    newSalinity[i, j] = Math.Max(0.0, Math.Min(salinityOcean, newSalinity[i, j]));
                    newTemperature[i, j] = Math.Max(0.0, Math.Min(temperatureOcean, newTemperature[i, j]));
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

            using (Brush waterBrush = new SolidBrush(Color.LightBlue))
            {
                g.FillRectangle(waterBrush, 0, 0, width, height);
            }

            double maxVorticity = 1.0;
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    int pixelX = i * width / gridPointsX;
                    int pixelY = height - (j + 1) * height / gridPointsZ;
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

            using (Pen velocityPen = new Pen(Color.Black, 1))
            {
                double maxVelocity = 0.5;
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
                        int endY = pixelY - arrowLengthY;
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