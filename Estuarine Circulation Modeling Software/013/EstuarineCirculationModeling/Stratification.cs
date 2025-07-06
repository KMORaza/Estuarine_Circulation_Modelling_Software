using System;
using System.Windows.Forms;
using System.Drawing;

namespace EstuarineCirculationModeling
{
    public class Stratification
    {
        private readonly int gridPoints;
        private readonly double dx; // Spatial step
        private readonly double estuaryDepth;
        private readonly double densityReference = 1000.0; // Reference density (kg/m³)
        private readonly double g = 9.81; // Gravitational acceleration (m/s²)
        private readonly BaroclinicFlow baroclinicFlow;
        private readonly PassiveScalarTransportEq passiveScalarTransport;
        private string turbulenceModel = "k-epsilon"; // Default turbulence model
        private double mixingParameter = 0.1; // Base mixing coefficient (m²/s)
        private double criticalRichardson = 0.25; // Critical Richardson number
        private double riverScalarConcentration = 1.0; // River boundary scalar concentration (kg/m³)
        private double[] richardsonNumber; // Gradient Richardson number at each grid point
        private double[] salinity; // Salinity profile for simulation
        private double[] temperature; // Temperature profile for simulation
        private double[] density; // Density profile for visualization
        private double[] passiveScalar; // Passive scalar concentration profile
        private Timer simulationTimer; // Timer for dynamic simulation
        private double simulationTime; // Current simulation time (s)
        private double dt = 3600.0; // Time step (1 hour)
        private readonly double salinityOcean = 35.0; // Ocean salinity (PSU)
        private readonly double temperatureOcean = 20.0; // Ocean temperature (°C)
        private readonly double scalarOcean = 0.0; // Ocean scalar concentration (kg/m³)
        private readonly double saltWedgePosition = 5000.0; // Salt wedge position (m)
        private readonly double estuaryLength = 10000.0; // Estuary length (m)

        public Stratification(double estuaryLength, double estuaryDepth, int gridPoints)
        {
            this.gridPoints = gridPoints;
            this.dx = estuaryLength / gridPoints;
            this.estuaryDepth = estuaryDepth;
            this.baroclinicFlow = new BaroclinicFlow(estuaryLength, estuaryDepth, gridPoints);
            this.passiveScalarTransport = new PassiveScalarTransportEq(estuaryLength, gridPoints);
            this.richardsonNumber = new double[gridPoints];
            this.salinity = new double[gridPoints];
            this.temperature = new double[gridPoints];
            this.density = new double[gridPoints];
            this.passiveScalar = new double[gridPoints];
            this.simulationTime = 0.0;

            // Initialize profiles
            ResetSimulation();
        }

        public void ComputeStratification(double[] velocity, double[] salinity, double[] temperature, double[] eddyViscosity)
        {
            // Calculate density profile
            for (int i = 0; i < gridPoints; i++)
            {
                density[i] = baroclinicFlow.CalculateDensity(salinity[i], temperature[i]);
                // Clamp density to prevent overflow
                density[i] = Math.Max(1000.0, Math.Min(1030.0, density[i]));
            }

            // Compute gradient Richardson number: Ri = (g/ρ)(∂ρ/∂z) / (∂u/∂z)²
            for (int i = 1; i < gridPoints - 1; i++)
            {
                double densityGradient = (density[i + 1] - density[i - 1]) / (2 * dx);
                double shear = Math.Pow((velocity[i + 1] - velocity[i - 1]) / (2 * dx), 2);
                richardsonNumber[i] = shear > 1e-10 ? (g / densityReference) * densityGradient / shear : 0.0;
                // Clamp Richardson number to prevent extreme values
                richardsonNumber[i] = Math.Max(-100.0, Math.Min(100.0, richardsonNumber[i]));
            }
            richardsonNumber[0] = richardsonNumber[1];
            richardsonNumber[gridPoints - 1] = richardsonNumber[gridPoints - 2];

            // Adjust mixing based on turbulence model
            switch (turbulenceModel)
            {
                case "k-epsilon":
                    for (int i = 0; i < gridPoints; i++)
                    {
                        eddyViscosity[i] *= (1.0 / (1.0 + Math.Max(0, richardsonNumber[i] / criticalRichardson)));
                    }
                    break;
                case "k-omega":
                    for (int i = 0; i < gridPoints; i++)
                    {
                        eddyViscosity[i] *= (1.0 + 0.5 * Math.Min(richardsonNumber[i], criticalRichardson)) / (1.0 + richardsonNumber[i]);
                    }
                    break;
                case "constant":
                    for (int i = 0; i < gridPoints; i++)
                    {
                        eddyViscosity[i] = mixingParameter;
                    }
                    break;
            }
            // Clamp eddy viscosity to prevent extreme values
            for (int i = 0; i < gridPoints; i++)
            {
                eddyViscosity[i] = Math.Max(0.0, Math.Min(1.0, eddyViscosity[i]));
            }
        }

        public double[] GetRichardsonNumber()
        {
            return (double[])richardsonNumber.Clone();
        }

        public double[] GetPassiveScalar()
        {
            return (double[])passiveScalar.Clone();
        }

        public void ShowStratificationWindow()
        {
            Form stratificationForm = new Form
            {
                Text = "Stratification Controls and Simulation",
                Size = new Size(600, 500),
                FormBorderStyle = FormBorderStyle.FixedDialog,
                MaximizeBox = false,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };

            // Control panel for inputs
            Panel controlPanel = new Panel
            {
                Location = new Point(10, 10),
                Size = new Size(250, 450),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };

            // Simulation panel
            Panel simulationPanel = new Panel
            {
                Location = new Point(270, 10),
                Size = new Size(300, 450),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = System.Drawing.Color.White
            };
            simulationPanel.Paint += (s, e) => DrawSimulation(s, e);

            // Turbulence model selection
            Label turbulenceLabel = new Label
            {
                Text = "Turbulence Model:",
                Location = new Point(10, 20),
                AutoSize = true,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };
            ComboBox turbulenceComboBox = new ComboBox
            {
                Location = new Point(10, 40),
                Size = new Size(200, 22),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.Color.White,
                ForeColor = System.Drawing.Color.Black,
                FlatStyle = FlatStyle.Flat
            };
            turbulenceComboBox.Items.AddRange(new string[] { "k-epsilon", "k-omega", "constant" });
            turbulenceComboBox.SelectedItem = turbulenceModel;
            turbulenceComboBox.SelectedIndexChanged += (s, e) => { turbulenceModel = turbulenceComboBox.SelectedItem.ToString(); UpdateSimulationParameters(); };

            // Mixing parameter
            Label mixingLabel = new Label
            {
                Text = "Mixing Coefficient (m²/s):",
                Location = new Point(10, 70),
                AutoSize = true,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };
            TextBox mixingTextBox = new TextBox
            {
                Location = new Point(10, 90),
                Size = new Size(200, 22),
                Text = mixingParameter.ToString("F2"),
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.Color.White,
                ForeColor = System.Drawing.Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            mixingTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(mixingTextBox.Text, out double value))
                {
                    mixingParameter = Math.Max(0.01, Math.Min(1.0, value));
                    UpdateSimulationParameters();
                }
            };

            // Critical Richardson number
            Label richardsonLabel = new Label
            {
                Text = "Critical Richardson Number:",
                Location = new Point(10, 120),
                AutoSize = true,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };
            TextBox richardsonTextBox = new TextBox
            {
                Location = new Point(10, 140),
                Size = new Size(200, 22),
                Text = criticalRichardson.ToString("F2"),
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.Color.White,
                ForeColor = System.Drawing.Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            richardsonTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(richardsonTextBox.Text, out double value))
                {
                    criticalRichardson = Math.Max(0.1, Math.Min(1.0, value));
                    UpdateSimulationParameters();
                }
            };

            // Passive scalar river concentration
            Label scalarLabel = new Label
            {
                Text = "River Scalar Conc. (kg/m³):",
                Location = new Point(10, 170),
                AutoSize = true,
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.SystemColors.Control,
                ForeColor = System.Drawing.Color.Black
            };
            TextBox scalarTextBox = new TextBox
            {
                Location = new Point(10, 190),
                Size = new Size(200, 22),
                Text = riverScalarConcentration.ToString("F2"),
                Font = new Font("Consolas", 9F),
                BackColor = System.Drawing.Color.White,
                ForeColor = System.Drawing.Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
            scalarTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(scalarTextBox.Text, out double value))
                {
                    riverScalarConcentration = Math.Max(0.0, Math.Min(10.0, value));
                    UpdateSimulationParameters();
                }
            };

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
                MessageBox.Show($"Turbulence Model: {turbulenceModel}\nMixing Coefficient: {mixingParameter:F2} m²/s\nCritical Ri: {criticalRichardson:F2}\nRiver Scalar: {riverScalarConcentration:F2} kg/m³", "Settings Applied", MessageBoxButtons.OK, MessageBoxIcon.Information);
                UpdateSimulationParameters();
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
                stratificationForm.Close();
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
                simulationTimer.Start();
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
                ResetSimulation();
                simulationPanel.Invalidate();
            };

            // Simulation timer
            simulationTimer = new Timer
            {
                Interval = 100 // Update every 100ms for smooth visualization
            };
            simulationTimer.Tick += (s, e) =>
            {
                UpdateSimulation();
                simulationPanel.Invalidate();
            };

            // Add controls to control panel
            controlPanel.Controls.AddRange(new Control[] { turbulenceLabel, turbulenceComboBox, mixingLabel, mixingTextBox, richardsonLabel, richardsonTextBox, scalarLabel, scalarTextBox, applyButton, closeButton, startButton, pauseButton, resetButton });

            // Add panels to form
            stratificationForm.Controls.Add(controlPanel);
            stratificationForm.Controls.Add(simulationPanel);

            // Start simulation
            simulationTimer.Start();

            // Show form
            stratificationForm.ShowDialog();
        }

        private void UpdateSimulationParameters()
        {
            // Reset simulation when parameters change
            ResetSimulation();
        }

        private void ResetSimulation()
        {
            // Initialize profiles (linear gradient from river to ocean)
            simulationTime = 0.0;
            for (int i = 0; i < gridPoints; i++)
            {
                double x = i * dx;
                double fraction = x / estuaryLength;
                salinity[i] = fraction * salinityOcean;
                temperature[i] = 15.0 + fraction * (temperatureOcean - 15.0);
                passiveScalar[i] = riverScalarConcentration * (1.0 - fraction); // Linear from river to ocean
                density[i] = baroclinicFlow.CalculateDensity(salinity[i], temperature[i]);
                // Clamp initial values
                salinity[i] = Math.Max(0.0, Math.Min(salinityOcean, salinity[i]));
                temperature[i] = Math.Max(0.0, Math.Min(temperatureOcean, temperature[i]));
                passiveScalar[i] = Math.Max(0.0, Math.Min(riverScalarConcentration, passiveScalar[i]));
                density[i] = Math.Max(1000.0, Math.Min(1030.0, density[i]));
            }
        }

        private void UpdateSimulation()
        {
            // Simplified velocity profile for simulation
            double[] velocity = new double[gridPoints];
            for (int i = 0; i < gridPoints; i++)
            {
                velocity[i] = 0.1 * (1.0 - i * dx / estuaryLength);
                // Clamp velocity to prevent extreme values
                velocity[i] = Math.Max(-1.0, Math.Min(1.0, velocity[i]));
            }

            // Update eddy viscosity based on current parameters
            double[] eddyViscosity = new double[gridPoints];
            for (int i = 0; i < gridPoints; i++)
            {
                eddyViscosity[i] = mixingParameter;
            }
            ComputeStratification(velocity, salinity, temperature, eddyViscosity);

            // Update salinity, temperature, and passive scalar
            salinity = baroclinicFlow.SolveSalinityTransport(salinity, velocity, dt, salinityOcean, saltWedgePosition, estuaryLength, eddyViscosity);
            temperature = baroclinicFlow.SolveTemperatureTransport(temperature, velocity, dt, temperatureOcean, saltWedgePosition, estuaryLength, eddyViscosity);
            passiveScalar = passiveScalarTransport.SolvePassiveScalarTransport(passiveScalar, velocity, dt, riverScalarConcentration, scalarOcean, saltWedgePosition, estuaryLength, eddyViscosity);

            // Clamp simulation outputs to prevent numerical instability
            for (int i = 0; i < gridPoints; i++)
            {
                salinity[i] = Math.Max(0.0, Math.Min(salinityOcean, salinity[i]));
                temperature[i] = Math.Max(0.0, Math.Min(temperatureOcean, temperature[i]));
                passiveScalar[i] = Math.Max(0.0, Math.Min(riverScalarConcentration, passiveScalar[i]));
                density[i] = Math.Max(1000.0, Math.Min(1030.0, density[i]));
            }

            // Advance simulation time
            simulationTime += dt;
        }

        private void DrawSimulation(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            Panel panel = (Panel)sender;
            g.Clear(Color.White);

            // Plot density and passive scalar profiles
            float panelWidth = panel.Width;
            float panelHeight = panel.Height;
            double maxDensity = 1030.0; // Max seawater density (kg/m³)
            double minDensity = 1000.0; // Min freshwater density (kg/m³)
            double maxScalar = Math.Max(1.0, riverScalarConcentration); // Adjust for user input
            double minScalar = 0.0; // Min scalar concentration

            // Prevent division by zero
            double densityRange = maxDensity - minDensity;
            double scalarRange = maxScalar - minScalar;
            if (densityRange < 1e-10) densityRange = 1e-10;
            if (scalarRange < 1e-10) scalarRange = 1e-10;

            // Draw axes
            g.DrawLine(Pens.Black, 30, 20, 30, panelHeight - 30); // Y-axis
            g.DrawLine(Pens.Black, 30, panelHeight - 30, panelWidth - 20, panelHeight - 30); // X-axis
            g.DrawString("Density (kg/m³) & Scalar (kg/m³)", new Font("Consolas", 9F), Brushes.Black, 10, 10);
            g.DrawString("Position (km)", new Font("Consolas", 9F), Brushes.Black, panelWidth - 50, panelHeight - 20);

            // Draw density profile (blue)
            PointF[] densityPoints = new PointF[gridPoints];
            for (int i = 0; i < gridPoints; i++)
            {
                float x = 30 + (i / (float)(gridPoints - 1)) * (panelWidth - 50);
                float normalizedDensity = (float)((density[i] - minDensity) / densityRange);
                normalizedDensity = Math.Max(0.0f, Math.Min(1.0f, normalizedDensity)); // Clamp to [0,1]
                float y = panelHeight - 30 - normalizedDensity * (panelHeight - 50);
                // Clamp y to panel bounds
                y = Math.Max(20, Math.Min(panelHeight - 30, y));
                densityPoints[i] = new PointF(x, y);
            }
            g.DrawLines(Pens.Blue, densityPoints);

            // Draw passive scalar profile (red)
            PointF[] scalarPoints = new PointF[gridPoints];
            for (int i = 0; i < gridPoints; i++)
            {
                float x = 30 + (i / (float)(gridPoints - 1)) * (panelWidth - 50);
                float normalizedScalar = (float)((passiveScalar[i] - minScalar) / scalarRange);
                normalizedScalar = Math.Max(0.0f, Math.Min(1.0f, normalizedScalar)); // Clamp to [0,1]
                float y = panelHeight - 30 - normalizedScalar * (panelHeight - 50);
                // Clamp y to panel bounds
                y = Math.Max(20, Math.Min(panelHeight - 30, y));
                scalarPoints[i] = new PointF(x, y);
            }
            g.DrawLines(Pens.Red, scalarPoints);

            // Draw labels for axes
            g.DrawString(maxDensity.ToString("F0"), new Font("Consolas", 9F), Brushes.Black, 5, 20);
            g.DrawString(minDensity.ToString("F0"), new Font("Consolas", 9F), Brushes.Black, 5, panelHeight - 30);
            g.DrawString("0", new Font("Consolas", 9F), Brushes.Black, 30, panelHeight - 20);
            g.DrawString((estuaryLength / 1000).ToString("F0"), new Font("Consolas", 9F), Brushes.Black, panelWidth - 30, panelHeight - 20);

            // Draw simulation time
            g.DrawString($"Time: {simulationTime / 3600:F1} hours", new Font("Consolas", 9F), Brushes.Black, panelWidth - 100, 10);

            // Draw legend
            g.DrawLine(Pens.Blue, panelWidth - 100, 30, panelWidth - 80, 30);
            g.DrawString("Density", new Font("Consolas", 9F), Brushes.Black, panelWidth - 70, 25);
            g.DrawLine(Pens.Red, panelWidth - 100, 50, panelWidth - 80, 50);
            g.DrawString("Scalar", new Font("Consolas", 9F), Brushes.Black, panelWidth - 70, 45);
        }
    }
}