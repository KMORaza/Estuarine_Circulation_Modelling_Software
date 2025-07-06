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
        private string turbulenceModel = "k-epsilon"; // Default turbulence model
        private double mixingParameter = 0.1; // Base mixing coefficient (m²/s)
        private double criticalRichardson = 0.25; // Critical Richardson number
        private double[] richardsonNumber; // Gradient Richardson number at each grid point
        private double[] salinity; // Salinity profile for simulation
        private double[] temperature; // Temperature profile for simulation
        private double[] density; // Density profile for visualization
        private Timer simulationTimer; // Timer for dynamic simulation
        private double simulationTime; // Current simulation time (s)
        private double dt = 3600.0; // Time step (1 hour)
        private readonly double salinityOcean = 35.0; // Ocean salinity (PSU)
        private readonly double temperatureOcean = 20.0; // Ocean temperature (°C)
        private readonly double saltWedgePosition = 5000.0; // Salt wedge position (m)
        private readonly double estuaryLength = 10000.0; // Estuary length (m)

        public Stratification(double estuaryLength, double estuaryDepth, int gridPoints)
        {
            this.gridPoints = gridPoints;
            this.dx = estuaryLength / gridPoints;
            this.estuaryDepth = estuaryDepth;
            this.baroclinicFlow = new BaroclinicFlow(estuaryLength, estuaryDepth, gridPoints);
            this.richardsonNumber = new double[gridPoints];
            this.salinity = new double[gridPoints];
            this.temperature = new double[gridPoints];
            this.density = new double[gridPoints];
            this.simulationTime = 0.0;

            // Initialize salinity and temperature profiles (linear gradient from river to ocean)
            ResetSimulation();
        }

        public void ComputeStratification(double[] velocity, double[] salinity, double[] temperature, double[] eddyViscosity)
        {
            // Calculate density profile
            for (int i = 0; i < gridPoints; i++)
            {
                density[i] = baroclinicFlow.CalculateDensity(salinity[i], temperature[i]);
            }

            // Compute gradient Richardson number: Ri = (g/ρ)(∂ρ/∂z) / (∂u/∂z)²
            for (int i = 1; i < gridPoints - 1; i++)
            {
                double densityGradient = (density[i + 1] - density[i - 1]) / (2 * dx); // Approximate vertical gradient
                double shear = Math.Pow((velocity[i + 1] - velocity[i - 1]) / (2 * dx), 2);
                richardsonNumber[i] = shear > 1e-10 ? (g / densityReference) * densityGradient / shear : 0.0;
            }
            richardsonNumber[0] = richardsonNumber[1];
            richardsonNumber[gridPoints - 1] = richardsonNumber[gridPoints - 2];

            // Adjust mixing based on turbulence model
            switch (turbulenceModel)
            {
                case "k-epsilon":
                    // k-ε model: Mixing suppressed when Ri > criticalRichardson
                    for (int i = 0; i < gridPoints; i++)
                    {
                        eddyViscosity[i] *= (1.0 / (1.0 + Math.Max(0, richardsonNumber[i] / criticalRichardson)));
                    }
                    break;
                case "k-omega":
                    // k-ω model: Enhanced mixing at low Ri, less sensitive to high Ri
                    for (int i = 0; i < gridPoints; i++)
                    {
                        eddyViscosity[i] *= (1.0 + 0.5 * Math.Min(richardsonNumber[i], criticalRichardson)) / (1.0 + richardsonNumber[i]);
                    }
                    break;
                case "constant":
                    // Constant diffusivity: No Ri dependence
                    for (int i = 0; i < gridPoints; i++)
                    {
                        eddyViscosity[i] = mixingParameter;
                    }
                    break;
            }
        }

        public double[] GetRichardsonNumber()
        {
            return (double[])richardsonNumber.Clone();
        }

        public void ShowStratificationWindow()
        {
            Form stratificationForm = new Form
            {
                Text = "Stratification Controls and Simulation",
                Size = new Size(600, 450),
                FormBorderStyle = FormBorderStyle.FixedSingle,
                MaximizeBox = false,
                Font = new Font("Tahoma", 8.25F)
            };

            // Control panel for inputs
            Panel controlPanel = new Panel
            {
                Location = new Point(10, 10),
                Size = new Size(250, 400),
                BorderStyle = BorderStyle.FixedSingle
            };

            // Simulation panel
            Panel simulationPanel = new Panel
            {
                Location = new Point(270, 10),
                Size = new Size(300, 400),
                BorderStyle = BorderStyle.FixedSingle
            };
            simulationPanel.Paint += (s, e) => DrawSimulation(s, e);

            // Turbulence model selection
            Label turbulenceLabel = new Label
            {
                Text = "Turbulence Model:",
                Location = new Point(10, 20),
                AutoSize = true
            };
            ComboBox turbulenceComboBox = new ComboBox
            {
                Location = new Point(10, 40),
                Size = new Size(200, 20),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            turbulenceComboBox.Items.AddRange(new string[] { "k-epsilon", "k-omega", "constant" });
            turbulenceComboBox.SelectedItem = turbulenceModel;
            turbulenceComboBox.SelectedIndexChanged += (s, e) => { turbulenceModel = turbulenceComboBox.SelectedItem.ToString(); UpdateSimulationParameters(); };

            // Mixing parameter
            Label mixingLabel = new Label
            {
                Text = "Mixing Coefficient (m²/s):",
                Location = new Point(10, 70),
                AutoSize = true
            };
            TextBox mixingTextBox = new TextBox
            {
                Location = new Point(10, 90),
                Size = new Size(200, 20),
                Text = mixingParameter.ToString("F2")
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
                AutoSize = true
            };
            TextBox richardsonTextBox = new TextBox
            {
                Location = new Point(10, 140),
                Size = new Size(200, 20),
                Text = criticalRichardson.ToString("F2")
            };
            richardsonTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(richardsonTextBox.Text, out double value))
                {
                    criticalRichardson = Math.Max(0.1, Math.Min(1.0, value));
                    UpdateSimulationParameters();
                }
            };

            // Apply button
            Button applyButton = new Button
            {
                Text = "Apply",
                Location = new Point(10, 180),
                Size = new Size(75, 23),
                Font = new Font("Tahoma", 8.25F)
            };
            applyButton.Click += (s, e) =>
            {
                MessageBox.Show($"Turbulence Model: {turbulenceModel}\nMixing Coefficient: {mixingParameter:F2} m²/s\nCritical Ri: {criticalRichardson:F2}", "Settings Applied", MessageBoxButtons.OK, MessageBoxIcon.Information);
                UpdateSimulationParameters();
            };

            // Close button
            Button closeButton = new Button
            {
                Text = "Close",
                Location = new Point(90, 180),
                Size = new Size(75, 23),
                Font = new Font("Tahoma", 8.25F)
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
                Location = new Point(10, 210),
                Size = new Size(75, 23),
                Font = new Font("Tahoma", 8.25F)
            };
            startButton.Click += (s, e) =>
            {
                simulationTimer.Start();
            };

            // Pause button
            Button pauseButton = new Button
            {
                Text = "Pause",
                Location = new Point(90, 210),
                Size = new Size(75, 23),
                Font = new Font("Tahoma", 8.25F)
            };
            pauseButton.Click += (s, e) =>
            {
                simulationTimer.Stop();
            };

            // Reset button
            Button resetButton = new Button
            {
                Text = "Reset",
                Location = new Point(170, 210),
                Size = new Size(75, 23),
                Font = new Font("Tahoma", 8.25F)
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
                simulationPanel.Invalidate(); // Redraw the simulation panel
            };

            // Add controls to control panel
            controlPanel.Controls.AddRange(new Control[] { turbulenceLabel, turbulenceComboBox, mixingLabel, mixingTextBox, richardsonLabel, richardsonTextBox, applyButton, closeButton, startButton, pauseButton, resetButton });

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
            // Initialize salinity and temperature profiles (linear gradient from river to ocean)
            simulationTime = 0.0;
            for (int i = 0; i < gridPoints; i++)
            {
                double x = i * dx;
                double fraction = x / estuaryLength;
                salinity[i] = fraction * salinityOcean;
                temperature[i] = 15.0 + fraction * (temperatureOcean - 15.0);
                density[i] = baroclinicFlow.CalculateDensity(salinity[i], temperature[i]);
            }
        }

        private void UpdateSimulation()
        {
            // Simplified velocity profile for simulation (assume steady flow for visualization)
            double[] velocity = new double[gridPoints];
            for (int i = 0; i < gridPoints; i++)
            {
                velocity[i] = 0.1 * (1.0 - i * dx / estuaryLength); // Linear decrease from river to ocean
            }

            // Update eddy viscosity based on current parameters
            double[] eddyViscosity = new double[gridPoints];
            for (int i = 0; i < gridPoints; i++)
            {
                eddyViscosity[i] = mixingParameter; // Initialize with base mixing
            }
            ComputeStratification(velocity, salinity, temperature, eddyViscosity);

            // Update salinity and temperature using BaroclinicFlow transport
            salinity = baroclinicFlow.SolveSalinityTransport(salinity, velocity, dt, salinityOcean, saltWedgePosition, estuaryLength, eddyViscosity);
            temperature = baroclinicFlow.SolveTemperatureTransport(temperature, velocity, dt, temperatureOcean, saltWedgePosition, estuaryLength, eddyViscosity);

            // Update density profile
            for (int i = 0; i < gridPoints; i++)
            {
                density[i] = baroclinicFlow.CalculateDensity(salinity[i], temperature[i]);
            }

            // Advance simulation time
            simulationTime += dt;
        }

        private void DrawSimulation(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            Panel panel = (Panel)sender;
            g.Clear(Color.White);

            // Plot density profile
            float panelWidth = panel.Width;
            float panelHeight = panel.Height;
            double maxDensity = 1030.0; // Typical max seawater density (kg/m³)
            double minDensity = 1000.0; // Typical min freshwater density (kg/m³)

            // Draw axes
            g.DrawLine(Pens.Black, 30, 20, 30, panelHeight - 30); // Y-axis
            g.DrawLine(Pens.Black, 30, panelHeight - 30, panelWidth - 20, panelHeight - 30); // X-axis
            g.DrawString("Density (kg/m³)", new Font("Tahoma", 8), Brushes.Black, 10, 10);
            g.DrawString("Position (km)", new Font("Tahoma", 8), Brushes.Black, panelWidth - 50, panelHeight - 20);

            // Draw density profile
            PointF[] points = new PointF[gridPoints];
            for (int i = 0; i < gridPoints; i++)
            {
                float x = 30 + (i / (float)(gridPoints - 1)) * (panelWidth - 50);
                float y = (float)(panelHeight - 30 - ((density[i] - minDensity) / (maxDensity - minDensity)) * (panelHeight - 50));
                points[i] = new PointF(x, y);
            }
            g.DrawLines(Pens.Blue, points);

            // Draw labels for axes
            g.DrawString(maxDensity.ToString("F0"), new Font("Tahoma", 8), Brushes.Black, 5, 20);
            g.DrawString(minDensity.ToString("F0"), new Font("Tahoma", 8), Brushes.Black, 5, panelHeight - 30);
            g.DrawString("0", new Font("Tahoma", 8), Brushes.Black, 30, panelHeight - 20);
            g.DrawString((estuaryLength / 1000).ToString("F0"), new Font("Tahoma", 8), Brushes.Black, panelWidth - 30, panelHeight - 20);

            // Draw simulation time
            g.DrawString($"Time: {simulationTime / 3600:F1} hours", new Font("Tahoma", 8), Brushes.Black, panelWidth - 100, 10);
        }
    }
}