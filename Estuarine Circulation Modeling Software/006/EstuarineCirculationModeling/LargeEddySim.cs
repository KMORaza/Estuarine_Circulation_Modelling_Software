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
                Size = new Size(600, 500),
                FormBorderStyle = FormBorderStyle.FixedDialog,
                MaximizeBox = false,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
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

        private static void SetupUI()
        {
            // Control panel
            Panel controlPanel = new Panel
            {
                Location = new Point(10, 10),
                Size = new Size(250, 450),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };

            // Simulation panel
            simulationPanel = new Panel
            {
                Location = new Point(270, 10),
                Size = new Size(300, 450),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = System.Drawing.Color.White
            };
            simulationPanel.Paint += (s, e) => DrawSimulation(e.Graphics, simulationPanel.Width, simulationPanel.Height);

            // Smagorinsky constant
            Label smagorinskyLabel = new Label
            {
                Text = "Smagorinsky Constant:",
                Location = new Point(10, 20),
                AutoSize = true,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };
            TextBox smagorinskyTextBox = new TextBox
            {
                Location = new Point(10, 40),
                Size = new Size(200, 22),
                Text = smagorinskyConstant.ToString("F2"),
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.Color.White,
                ForeColor = System.Drawing.Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            smagorinskyTextBox.TextChanged += (s, e) => ValidateSmagorinsky(smagorinskyTextBox);

            // Grid resolution
            Label gridLabel = new Label
            {
                Text = "Grid Points (XxZ):",
                Location = new Point(10, 70),
                AutoSize = true,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };
            TextBox gridTextBox = new TextBox
            {
                Location = new Point(10, 90),
                Size = new Size(200, 22),
                Text = $"{gridPointsX}x{gridPointsZ}",
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.Color.White,
                ForeColor = System.Drawing.Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            gridTextBox.TextChanged += (s, e) => ValidateGrid(gridTextBox);

            // River inflow
            Label inflowLabel = new Label
            {
                Text = "River Inflow (m³/s):",
                Location = new Point(10, 120),
                AutoSize = true,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };
            TextBox inflowTextBox = new TextBox
            {
                Location = new Point(10, 140),
                Size = new Size(200, 22),
                Text = riverInflow.ToString("F2"),
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.Color.White,
                ForeColor = System.Drawing.Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            inflowTextBox.TextChanged += (s, e) => ValidateInflow(inflowTextBox);

            // Tidal amplitude
            Label tidalLabel = new Label
            {
                Text = "Tidal Amplitude (m):",
                Location = new Point(10, 170),
                AutoSize = true,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };
            TextBox tidalTextBox = new TextBox
            {
                Location = new Point(10, 190),
                Size = new Size(200, 22),
                Text = tidalAmplitude.ToString("F2"),
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.Color.White,
                ForeColor = System.Drawing.Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            tidalTextBox.TextChanged += (s, e) => ValidateTidal(tidalTextBox);

            // Apply button
            Button applyButton = new Button
            {
                Text = "Apply",
                Location = new Point(10, 230),
                Size = new Size(75, 23),
                Font = new Font("Consolas", 9F),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
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
                Location = new Point(90, 230),
                Size = new Size(75, 23),
                Font = new Font("Consolas", 9F),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
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
                Location = new Point(10, 260),
                Size = new Size(75, 23),
                Font = new Font("Consolas", 9F),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };
            startButton.Click += (s, e) =>
            {
                if (ValidateAllInputs(smagorinskyTextBox, gridTextBox, inflowTextBox, tidalTextBox))
                {
                    simulationTimer.Start();
                    isSimulationRunning = true;
                    UpdateStatus("Running");
                }
            };

            // Pause button
            Button pauseButton = new Button
            {
                Text = "Pause",
                Location = new Point(90, 260),
                Size = new Size(75, 23),
                Font = new Font("Consolas", 9F),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
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
                Location = new Point(170, 260),
                Size = new Size(75, 23),
                Font = new Font("Consolas", 9F),
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
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
                Location = new Point(10, 290),
                Size = new Size(200, 40),
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };

            // Simulation timer
            simulationTimer = new Timer
            {
                Interval = 100 // Update every 100ms
            };
            simulationTimer.Tick += (s, e) =>
            {
                UpdateSimulation();
                simulationPanel.Invalidate();
                UpdateStatus($"Running (Time: {simulationTime:F1} s)");
            };

            // Add controls to control panel
            controlPanel.Controls.AddRange(new Control[] { smagorinskyLabel, smagorinskyTextBox, gridLabel, gridTextBox, inflowLabel, inflowTextBox, tidalLabel, tidalTextBox, applyButton, closeButton, startButton, pauseButton, resetButton, statusLabel });

            // Add panels to form
            lesForm.Controls.Add(controlPanel);
            lesForm.Controls.Add(simulationPanel);
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
            if (textBox.Text.Contains("x"))
            {
                var parts = textBox.Text.Split('x');
                if (parts.Length == 2 && int.TryParse(parts[0], out int x) && int.TryParse(parts[1], out int z))
                {
                    gridPointsX = Math.Max(20, Math.Min(100, x));
                    gridPointsZ = Math.Max(10, Math.Min(50, z));
                    textBox.Text = $"{gridPointsX}x{gridPointsZ}";
                    UpdateStatus("Ready");
                    return true;
                }
            }
            UpdateStatus("Invalid Grid Format (e.g., 50x20)");
            return false;
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
            for (int i = 1; i < gridPointsX - 1; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    pressureGradientX[i, j] = g * density[i, j] * (density[i + 1, j] - density[i - 1, j]) / (2 * dx * densityReference);
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

            // Solve momentum equations (simplified Navier-Stokes with LES)
            for (int i = 1; i < gridPointsX - 1; i++)
            {
                for (int j = 1; j < gridPointsZ - 1; j++)
                {
                    double advectionU = velocityU[i, j] * (velocityU[i + 1, j] - velocityU[i - 1, j]) / (2 * dx) +
                                       velocityW[i, j] * (velocityU[i, j + 1] - velocityU[i, j - 1]) / (2 * dz);
                    double diffusionU = eddyViscosity[i, j] * ((velocityU[i + 1, j] - 2 * velocityU[i, j] + velocityU[i - 1, j]) / (dx * dx) +
                                                              (velocityU[i, j + 1] - 2 * velocityU[i, j] + velocityU[i, j - 1]) / (dz * dz));
                    newVelocityU[i, j] = velocityU[i, j] + dt * (-advectionU - pressureGradientX[i, j] / densityReference + diffusionU);

                    double advectionW = velocityU[i, j] * (velocityW[i + 1, j] - velocityW[i - 1, j]) / (2 * dx) +
                                       velocityW[i, j] * (velocityW[i, j + 1] - velocityW[i, j - 1]) / (2 * dz);
                    double diffusionW = eddyViscosity[i, j] * ((velocityW[i + 1, j] - 2 * velocityW[i, j] + velocityW[i - 1, j]) / (dx * dx) +
                                                              (velocityW[i, j + 1] - 2 * velocityW[i, j] + velocityW[i, j - 1]) / (dz * dz));
                    newVelocityW[i, j] = velocityW[i, j] + dt * (-advectionW + diffusionW);

                    // Clamp velocities
                    newVelocityU[i, j] = Math.Max(-1.0, Math.Min(1.0, newVelocityU[i, j]));
                    newVelocityW[i, j] = Math.Max(-0.1, Math.Min(0.1, newVelocityW[i, j]));
                }
            }

            // Boundary conditions
            for (int j = 0; j < gridPointsZ; j++)
            {
                newVelocityU[0, j] = riverInflow / estuaryDepth;
                newVelocityU[gridPointsX - 1, j] = tidalAmplitude * Math.Cos(2 * Math.PI * simulationTime / tidalPeriod);
                newVelocityW[0, j] = 0.0;
                newVelocityW[gridPointsX - 1, j] = 0.0;
                newSalinity[0, j] = 0.0;
                newSalinity[gridPointsX - 1, j] = salinityOcean;
                newTemperature[0, j] = 15.0;
                newTemperature[gridPointsX - 1, j] = temperatureOcean;
            }
            for (int i = 0; i < gridPointsX; i++)
            {
                newVelocityW[i, 0] = 0.0; // Bottom boundary
                newVelocityW[i, gridPointsZ - 1] = 0.0; // Surface boundary
            }

            // Update salinity and temperature (advection-diffusion)
            for (int i = 1; i < gridPointsX - 1; i++)
            {
                for (int j = 1; j < gridPointsZ - 1; j++)
                {
                    double advectionS = velocityU[i, j] * (salinity[i + 1, j] - salinity[i - 1, j]) / (2 * dx) +
                                       velocityW[i, j] * (salinity[i, j + 1] - salinity[i, j - 1]) / (2 * dz);
                    double diffusionS = eddyViscosity[i, j] * ((salinity[i + 1, j] - 2 * salinity[i, j] + salinity[i - 1, j]) / (dx * dx) +
                                                              (salinity[i, j + 1] - 2 * salinity[i, j] + salinity[i, j - 1]) / (dz * dz));
                    newSalinity[i, j] = salinity[i, j] + dt * (-advectionS + diffusionS);

                    double advectionT = velocityU[i, j] * (temperature[i + 1, j] - temperature[i - 1, j]) / (2 * dx) +
                                       velocityW[i, j] * (temperature[i, j + 1] - temperature[i, j - 1]) / (2 * dz);
                    double diffusionT = eddyViscosity[i, j] * ((temperature[i + 1, j] - 2 * temperature[i, j] + temperature[i - 1, j]) / (dx * dx) +
                                                              (temperature[i, j + 1] - 2 * temperature[i, j] + temperature[i, j - 1]) / (dz * dz));
                    newTemperature[i, j] = temperature[i, j] + dt * (-advectionT + diffusionT);

                    // Clamp scalars
                    newSalinity[i, j] = Math.Max(0.0, Math.Min(salinityOcean, newSalinity[i, j]));
                    newTemperature[i, j] = Math.Max(0.0, Math.Min(temperatureOcean, newTemperature[i, j]));
                }
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

            // Draw labels
            using (Font font = new Font("Consolas", 9F))
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