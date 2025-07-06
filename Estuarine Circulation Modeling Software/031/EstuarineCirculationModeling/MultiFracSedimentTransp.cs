using System;
using System.Collections.Generic;
using System.Windows.Forms;
using System.Drawing;

namespace EstuarineCirculationModeling
{
    public class MultiFracSedimentTransp : Form
    {
        private TabControl tabControl;
        private TabPage visualizationTab;
        private TabPage hovmoellerTab;
        private TabPage budgetMapTab;
        private Panel visualizationPanel;
        private Panel hovmoellerPanel;
        private Panel budgetMapPanel;
        private TextBox sandConcentrationTextBox;
        private TextBox siltConcentrationTextBox;
        private TextBox clayConcentrationTextBox;
        private TextBox tidalAmplitudeTextBox;
        private TextBox tidalPeriodTextBox;
        private TextBox riverInflowTextBox;
        private ComboBox sedimentSelector;
        private CheckBox showHovmoellerCheckBox;
        private CheckBox showBudgetMapCheckBox;
        private Label settlingVelocitiesLabel;
        private Button startSimButton;
        private Button pauseSimButton;
        private Button resetSimButton;
        private TextBox outputTextBox;
        private Timer simulationTimer;
        private bool isRunning = false;
        private double[][] sedimentConcentrations; // [fraction][gridPoint]
        private readonly string[] sedimentTypes = { "Sand", "Silt", "Clay" };
        private readonly double[] grainSizes = { 0.0002, 0.00002, 0.000002 }; // m (200 µm, 20 µm, 2 µm)
        private readonly double[] settlingVelocities; // m/s
        private readonly double[] diffusionCoefficients = { 0.001, 0.0005, 0.0002 }; // m²/s
        private readonly double[] criticalShearStresses = { 0.2, 0.1, 0.05 }; // N/m² (for erosion)
        private double[] positions;
        private double[] velocityField; // m/s
        private List<double[]> sedimentFluxes; // [time][fraction * gridPoint] (net flux: deposition - erosion)
        private double time = 0.0;
        private double estuaryLength = 10000.0; // m
        private double estuaryDepth = 10.0; // m
        private double timeStep = 0.1; // s
        private int gridPoints = 100;
        private double tidalAmplitude = 1.0; // m
        private double tidalPeriod = 600.0; // s (~10 min)
        private double riverInflow = 0.1; // m³/s
        private const int maxTimeSteps = 6000; // 600 s / 0.1 s

        public MultiFracSedimentTransp()
        {
            settlingVelocities = CalculateSettlingVelocities();
            InitializeComponents();
            InitializeSimulation();
        }

        private double[] CalculateSettlingVelocities()
        {
            double[] velocities = new double[grainSizes.Length];
            double gravity = 9.81; // m/s²
            double waterViscosity = 1.0e-6; // m²/s (kinematic viscosity at 20°C)
            double sedimentDensity = 2650.0; // kg/m³ (quartz)
            double waterDensity = 1000.0; // kg/m³
            double reductionFactor = 0.001; // Reduced for sand transport

            for (int i = 0; i < grainSizes.Length; i++)
            {
                velocities[i] = reductionFactor * (sedimentDensity - waterDensity) * gravity * grainSizes[i] * grainSizes[i] / (18 * waterViscosity);
            }
            return velocities;
        }

        private void InitializeComponents()
        {
            this.Text = "Multi-Fraction Sediment Transport Simulation";
            this.Size = new Size(1200, 900);
            this.FormBorderStyle = FormBorderStyle.Sizable;
            this.MinimumSize = new Size(1000, 700);
            this.Font = new Font("Verdana", 10F);

            var controlPanel = new Panel
            {
                Location = new Point(10, 10),
                Size = new Size(300, 840),
                BorderStyle = BorderStyle.FixedSingle,
                Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Bottom
            };

            // Sediment Inputs Group
            var sedimentGroupLabel = new Label
            {
                Location = new Point(10, 10),
                Size = new Size(280, 20),
                Text = "Sediment Concentrations",
                Font = new Font("Verdana", 10F, FontStyle.Bold)
            };
            var sandLabel = new Label
            {
                Location = new Point(10, 40),
                Size = new Size(280, 20),
                Text = "Sand (kg/m³):"
            };
            sandConcentrationTextBox = new TextBox
            {
                Location = new Point(10, 60),
                Size = new Size(260, 25),
                Text = "0.1"
            };
            sandConcentrationTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(sandConcentrationTextBox, "Initial sand concentration (kg/m³)");
            var siltLabel = new Label
            {
                Location = new Point(10, 90),
                Size = new Size(280, 20),
                Text = "Silt (kg/m³):"
            };
            siltConcentrationTextBox = new TextBox
            {
                Location = new Point(10, 110),
                Size = new Size(260, 25),
                Text = "0.05"
            };
            siltConcentrationTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(siltConcentrationTextBox, "Initial silt concentration (kg/m³)");
            var clayLabel = new Label
            {
                Location = new Point(10, 140),
                Size = new Size(280, 20),
                Text = "Clay (kg/m³):"
            };
            clayConcentrationTextBox = new TextBox
            {
                Location = new Point(10, 160),
                Size = new Size(260, 25),
                Text = "0.02"
            };
            clayConcentrationTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(clayConcentrationTextBox, "Initial clay concentration (kg/m³)");

            // Flow Inputs Group
            var flowGroupLabel = new Label
            {
                Location = new Point(10, 200),
                Size = new Size(280, 20),
                Text = "Flow Parameters",
                Font = new Font("Verdana", 10F, FontStyle.Bold)
            };
            var tidalAmplitudeLabel = new Label
            {
                Location = new Point(10, 230),
                Size = new Size(280, 20),
                Text = "Tidal Amplitude (m):"
            };
            tidalAmplitudeTextBox = new TextBox
            {
                Location = new Point(10, 250),
                Size = new Size(260, 25),
                Text = "1.0"
            };
            tidalAmplitudeTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(tidalAmplitudeTextBox, "Tidal wave amplitude (0-10 m)");
            var tidalPeriodLabel = new Label
            {
                Location = new Point(10, 280),
                Size = new Size(280, 20),
                Text = "Tidal Period (s):"
            };
            tidalPeriodTextBox = new TextBox
            {
                Location = new Point(10, 300),
                Size = new Size(260, 25),
                Text = "600"
            };
            tidalPeriodTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(tidalPeriodTextBox, "Tidal period (600-86400 s)");
            var riverInflowLabel = new Label
            {
                Location = new Point(10, 330),
                Size = new Size(280, 20),
                Text = "River Inflow (m³/s):"
            };
            riverInflowTextBox = new TextBox
            {
                Location = new Point(10, 350),
                Size = new Size(260, 25),
                Text = "0.1"
            };
            riverInflowTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(riverInflowTextBox, "River inflow rate (0-100 m³/s)");

            // Visualization Options Group
            var vizGroupLabel = new Label
            {
                Location = new Point(10, 390),
                Size = new Size(280, 20),
                Text = "Visualization Options",
                Font = new Font("Verdana", 10F, FontStyle.Bold)
            };
            var sedimentSelectorLabel = new Label
            {
                Location = new Point(10, 420),
                Size = new Size(280, 20),
                Text = "Sediment for Maps:"
            };
            sedimentSelector = new ComboBox
            {
                Location = new Point(10, 440),
                Size = new Size(260, 25),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            sedimentSelector.Items.AddRange(new object[] { "Sand", "Silt", "Clay", "All" });
            sedimentSelector.SelectedIndex = 0;
            sedimentSelector.SelectedIndexChanged += (s, e) => { hovmoellerPanel.Invalidate(); budgetMapPanel.Invalidate(); };
            sedimentSelector.MouseHover += (s, e) => new ToolTip().SetToolTip(sedimentSelector, "Select sediment type for Hovmöller and budget maps");

            showHovmoellerCheckBox = new CheckBox
            {
                Location = new Point(10, 470),
                Size = new Size(280, 25),
                Text = "Show Hovmöller Diagram",
                Checked = true
            };
            showHovmoellerCheckBox.MouseHover += (s, e) => new ToolTip().SetToolTip(showHovmoellerCheckBox, "Toggle temporal sediment flux diagram");
            showBudgetMapCheckBox = new CheckBox
            {
                Location = new Point(10, 500),
                Size = new Size(280, 25),
                Text = "Show Sediment Budget Map",
                Checked = true
            };
            showBudgetMapCheckBox.MouseHover += (s, e) => new ToolTip().SetToolTip(showBudgetMapCheckBox, "Toggle spatial sediment budget map");

            settlingVelocitiesLabel = new Label
            {
                Location = new Point(10, 530),
                Size = new Size(280, 80),
                Text = $"Settling Velocities (m/s):\nSand: {settlingVelocities[0]:F6}\nSilt: {settlingVelocities[1]:F6}\nClay: {settlingVelocities[2]:F6}"
            };

            // Control Buttons
            startSimButton = new Button
            {
                Location = new Point(10, 620),
                Size = new Size(120, 30),
                Text = "Start/Resume",
                FlatStyle = FlatStyle.Flat
            };
            startSimButton.Click += StartSimButton_Click;
            startSimButton.MouseHover += (s, e) => new ToolTip().SetToolTip(startSimButton, "Start or resume the simulation");
            pauseSimButton = new Button
            {
                Location = new Point(140, 620),
                Size = new Size(120, 30),
                Text = "Pause",
                FlatStyle = FlatStyle.Flat,
                Enabled = false
            };
            pauseSimButton.Click += PauseSimButton_Click;
            pauseSimButton.MouseHover += (s, e) => new ToolTip().SetToolTip(pauseSimButton, "Pause the simulation");
            resetSimButton = new Button
            {
                Location = new Point(10, 660),
                Size = new Size(120, 30),
                Text = "Reset",
                FlatStyle = FlatStyle.Flat
            };
            resetSimButton.Click += ResetSimButton_Click;
            resetSimButton.MouseHover += (s, e) => new ToolTip().SetToolTip(resetSimButton, "Reset the simulation to initial conditions");

            outputTextBox = new TextBox
            {
                Location = new Point(320, 710),
                Size = new Size(860, 150),
                Multiline = true,
                ReadOnly = true,
                ScrollBars = ScrollBars.Vertical,
                Font = new Font("Courier New", 10F),
                BackColor = Color.Black,
                ForeColor = Color.White,
                Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Right | AnchorStyles.Bottom
            };

            tabControl = new TabControl
            {
                Location = new Point(320, 10),
                Size = new Size(860, 690),
                Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Right | AnchorStyles.Bottom
            };

            visualizationTab = new TabPage { Text = "Profiles" };
            visualizationPanel = new Panel
            {
                Location = new Point(0, 0),
                Size = new Size(856, 666),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            visualizationPanel.Paint += VisualizationPanel_Paint;
            visualizationTab.Controls.Add(visualizationPanel);

            hovmoellerTab = new TabPage { Text = "Hovmöller" };
            hovmoellerPanel = new Panel
            {
                Location = new Point(0, 0),
                Size = new Size(856, 666),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White,
                Visible = showHovmoellerCheckBox.Checked
            };
            hovmoellerPanel.Paint += HovmoellerPanel_Paint;
            showHovmoellerCheckBox.CheckedChanged += (s, e) => { hovmoellerTab.Enabled = showHovmoellerCheckBox.Checked; if (!showHovmoellerCheckBox.Checked && tabControl.SelectedTab == hovmoellerTab) tabControl.SelectedTab = visualizationTab; };
            hovmoellerTab.Controls.Add(hovmoellerPanel);

            budgetMapTab = new TabPage { Text = "Budget Map" };
            budgetMapPanel = new Panel
            {
                Location = new Point(0, 0),
                Size = new Size(856, 666),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White,
                Visible = showBudgetMapCheckBox.Checked
            };
            budgetMapPanel.Paint += BudgetMapPanel_Paint;
            showBudgetMapCheckBox.CheckedChanged += (s, e) => { budgetMapTab.Enabled = showBudgetMapCheckBox.Checked; if (!showBudgetMapCheckBox.Checked && tabControl.SelectedTab == budgetMapTab) tabControl.SelectedTab = visualizationTab; };
            budgetMapTab.Controls.Add(budgetMapPanel);

            tabControl.TabPages.AddRange(new[] { visualizationTab, hovmoellerTab, budgetMapTab });

            controlPanel.Controls.AddRange(new Control[] {
                sedimentGroupLabel, sandLabel, sandConcentrationTextBox,
                siltLabel, siltConcentrationTextBox,
                clayLabel, clayConcentrationTextBox,
                flowGroupLabel, tidalAmplitudeLabel, tidalAmplitudeTextBox,
                tidalPeriodLabel, tidalPeriodTextBox,
                riverInflowLabel, riverInflowTextBox,
                vizGroupLabel, sedimentSelectorLabel, sedimentSelector,
                showHovmoellerCheckBox, showBudgetMapCheckBox,
                settlingVelocitiesLabel,
                startSimButton, pauseSimButton, resetSimButton
            });

            this.Controls.AddRange(new Control[] { controlPanel, outputTextBox, tabControl });
        }

        private void InitializeSimulation()
        {
            double[] initialConcentrations = ParseInputs();
            sedimentConcentrations = new double[sedimentTypes.Length][];
            for (int f = 0; f < sedimentTypes.Length; f++)
            {
                sedimentConcentrations[f] = new double[gridPoints];
                for (int i = 0; i < gridPoints; i++)
                    sedimentConcentrations[f][i] = initialConcentrations[f];
            }
            positions = new double[gridPoints];
            velocityField = new double[gridPoints];
            sedimentFluxes = new List<double[]>();
            for (int i = 0; i < gridPoints; i++)
            {
                positions[i] = i * estuaryLength / (gridPoints - 1);
                velocityField[i] = riverInflow / estuaryDepth;
            }

            simulationTimer = new Timer
            {
                Interval = 100
            };
            simulationTimer.Tick += (s, e) => UpdateSimulation();
            UpdateButtonStates();
        }

        private double[] ParseInputs()
        {
            try
            {
                double sand = double.Parse(sandConcentrationTextBox.Text.Trim());
                double silt = double.Parse(siltConcentrationTextBox.Text.Trim());
                double clay = double.Parse(clayConcentrationTextBox.Text.Trim());
                tidalAmplitude = double.Parse(tidalAmplitudeTextBox.Text.Trim());
                tidalPeriod = double.Parse(tidalPeriodTextBox.Text.Trim());
                riverInflow = double.Parse(riverInflowTextBox.Text.Trim());

                if (sand < 0 || silt < 0 || clay < 0)
                    throw new Exception("Concentrations must be non-negative.");
                if (tidalAmplitude < 0 || tidalAmplitude > 10.0)
                    throw new Exception("Tidal amplitude must be between 0 and 10 m.");
                if (tidalPeriod < 600.0 || tidalPeriod > 86400.0)
                    throw new Exception("Tidal period must be between 600 and 86400 s.");
                if (riverInflow < 0 || riverInflow > 100.0)
                    throw new Exception("River inflow must be between 0 and 100 m³/s.");

                return new double[] { sand, silt, clay };
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Invalid input: {ex.Message}", "Input Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                return new double[] { 0.1, 0.05, 0.02 };
            }
        }

        private void StartSimButton_Click(object sender, EventArgs e)
        {
            try
            {
                if (!isRunning)
                {
                    ParseInputs();
                    UpdateVelocityField();
                    isRunning = true;
                    simulationTimer.Start();
                    outputTextBox.AppendText(sedimentFluxes.Count == 0 ? $"Simulation started with {sedimentTypes.Length} sediment types.\r\n" : "Simulation resumed.\r\n");
                    UpdateButtonStates();
                    visualizationPanel.Invalidate();
                    hovmoellerPanel.Invalidate();
                    budgetMapPanel.Invalidate();
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Error: {ex.Message}", "Input Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private void PauseSimButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isRunning = false;
            outputTextBox.AppendText("Simulation paused.\r\n");
            UpdateButtonStates();
        }

        private void ResetSimButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isRunning = false;
            time = 0.0;
            double[] initialConcentrations = ParseInputs();
            sedimentConcentrations = new double[sedimentTypes.Length][];
            for (int f = 0; f < sedimentTypes.Length; f++)
            {
                sedimentConcentrations[f] = new double[gridPoints];
                for (int i = 0; i < gridPoints; i++)
                    sedimentConcentrations[f][i] = initialConcentrations[f];
            }
            sedimentFluxes.Clear();
            UpdateVelocityField();
            outputTextBox.Clear();
            outputTextBox.AppendText("Simulation reset.\r\n");
            visualizationPanel.Invalidate();
            hovmoellerPanel.Invalidate();
            budgetMapPanel.Invalidate();
            UpdateButtonStates();
        }

        private void UpdateVelocityField()
        {
            double dx = estuaryLength / (gridPoints - 1);
            double tidalVelocity = 2.0 * tidalAmplitude * (2 * Math.PI / tidalPeriod) * Math.Sin(2 * Math.PI * time / tidalPeriod);
            double riverVelocity = riverInflow / estuaryDepth;
            for (int i = 0; i < gridPoints; i++)
            {
                double x = positions[i];
                double distanceFactor = 1.0 - x / estuaryLength;
                velocityField[i] = riverVelocity + tidalVelocity * distanceFactor;
            }
        }

        private void UpdateSimulation()
        {
            UpdateVelocityField();
            double dx = estuaryLength / (gridPoints - 1);
            double[] currentFluxes = new double[sedimentTypes.Length * gridPoints];
            double waterDensity = 1000.0; // kg/m³
            double dragCoefficient = 0.0025; // Dimensionless

            for (int f = 0; f < sedimentTypes.Length; f++)
            {
                double[] newConcentrations = new double[gridPoints];
                for (int i = 1; i < gridPoints - 1; i++)
                {
                    double advection = velocityField[i] >= 0
                        ? -velocityField[i] * (sedimentConcentrations[f][i] - sedimentConcentrations[f][i - 1]) / dx
                        : -velocityField[i] * (sedimentConcentrations[f][i + 1] - sedimentConcentrations[f][i]) / dx;
                    double diffusion = diffusionCoefficients[f] * (sedimentConcentrations[f][i + 1] - 2 * sedimentConcentrations[f][i] + sedimentConcentrations[f][i - 1]) / (dx * dx);
                    double settling = -settlingVelocities[f] * sedimentConcentrations[f][i];
                    newConcentrations[i] = sedimentConcentrations[f][i] + timeStep * (advection + diffusion + settling);
                    if (newConcentrations[i] < 0) newConcentrations[i] = 0;

                    double depositionRate = settlingVelocities[f] * sedimentConcentrations[f][i];
                    double shearStress = waterDensity * dragCoefficient * velocityField[i] * velocityField[i];
                    double erosionRate = shearStress > criticalShearStresses[f] ? 0.001 * (shearStress - criticalShearStresses[f]) : 0.0;
                    currentFluxes[f * gridPoints + i] = depositionRate - erosionRate;
                }
                newConcentrations[0] = sedimentConcentrations[f][0];
                newConcentrations[gridPoints - 1] = 0.0;
                sedimentConcentrations[f] = newConcentrations;
            }
            sedimentFluxes.Add(currentFluxes);
            if (sedimentFluxes.Count > maxTimeSteps)
                sedimentFluxes.RemoveAt(0);

            time += timeStep;
            visualizationPanel.Invalidate();
            hovmoellerPanel.Invalidate();
            budgetMapPanel.Invalidate();
            UpdateOutput();
        }

        private void UpdateOutput()
        {
            double tidalPhase = (2 * Math.PI * time / tidalPeriod) % (2 * Math.PI);
            string phaseStr = tidalPhase < Math.PI ? "Flood" : "Ebb";
            double maxVelocity = velocityField[0];
            double minVelocity = velocityField[gridPoints - 1];
            string output = $"Time: {time:F2}s | Tidal Phase: {phaseStr} | Vel (mid): {velocityField[gridPoints / 2]:F6} m/s | Vel (range): [{minVelocity:F6}, {maxVelocity:F6}] m/s | Conc (mid): ";
            for (int f = 0; f < sedimentTypes.Length; f++)
            {
                double flux = velocityField[gridPoints / 2] * sedimentConcentrations[f][gridPoints / 2];
                output += $"{sedimentTypes[f]}: {sedimentConcentrations[f][gridPoints / 2]:F6} kg/m³ (Flux: {flux:F6} kg/m²s) ";
            }
            outputTextBox.AppendText(output + "\r\n");
        }

        private void VisualizationPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            g.Clear(Color.White);

            int panelWidth = visualizationPanel.Width;
            int panelHeight = visualizationPanel.Height;
            Color[] colors = { Color.Red, Color.Blue, Color.Green };

            for (int f = 0; f < sedimentTypes.Length; f++)
            {
                using (Pen pen = new Pen(colors[f], 2))
                {
                    for (int i = 0; i < gridPoints - 1; i++)
                    {
                        float x1 = i * panelWidth / (float)(gridPoints - 1);
                        float x2 = (i + 1) * panelWidth / (float)(gridPoints - 1);
                        float y1 = panelHeight - (float)(sedimentConcentrations[f][i] * panelHeight / 0.2);
                        float y2 = panelHeight - (float)(sedimentConcentrations[f][i + 1] * panelHeight / 0.2);
                        g.DrawLine(pen, x1, y1, x2, y2);
                    }
                }
            }

            using (Pen velocityPen = new Pen(Color.Black, 2))
            {
                double maxVelocity = 0.2;
                for (int i = 0; i < gridPoints - 1; i++)
                {
                    float x1 = i * panelWidth / (float)(gridPoints - 1);
                    float x2 = (i + 1) * panelWidth / (float)(gridPoints - 1);
                    float y1 = panelHeight / 2 - (float)(velocityField[i] * panelHeight / (2 * maxVelocity));
                    float y2 = panelHeight / 2 - (float)(velocityField[i + 1] * panelHeight / (2 * maxVelocity));
                    g.DrawLine(velocityPen, x1, y1, x2, y2);
                }
            }

            g.DrawLine(Pens.Black, 0, panelHeight - 10, panelWidth, panelHeight - 10);
            g.DrawLine(Pens.Black, 10, 0, 10, panelHeight);
            g.DrawString("Distance (m)", new Font("Verdana", 10), Brushes.Black, panelWidth - 120, panelHeight - 30);
            g.DrawString("Concentration (kg/m³)", new Font("Verdana", 10), Brushes.Black, 10, 10);
            g.DrawString("Velocity (m/s)", new Font("Verdana", 10), Brushes.Black, 10, panelHeight / 2 - 20);

            for (int f = 0; f < sedimentTypes.Length; f++)
            {
                g.FillRectangle(new SolidBrush(colors[f]), 20, 20 + f * 25, 15, 15);
                g.DrawString(sedimentTypes[f], new Font("Verdana", 10), Brushes.Black, 40, 20 + f * 25);
            }
            g.FillRectangle(new SolidBrush(Color.Black), 20, 95, 15, 15);
            g.DrawString("Velocity", new Font("Verdana", 10), Brushes.Black, 40, 95);

            double tidalPhase = (2 * Math.PI * time / tidalPeriod) % (2 * Math.PI);
            string phaseStr = tidalPhase < Math.PI ? "Flood" : "Ebb";
            g.DrawString($"Time: {time:F2} s | Phase: {phaseStr}", new Font("Verdana", 10), Brushes.Black, panelWidth - 180, 10);
        }

        private void HovmoellerPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            g.Clear(Color.White);

            int panelWidth = hovmoellerPanel.Width;
            int panelHeight = hovmoellerPanel.Height;
            int selectedSediment = sedimentSelector.SelectedIndex == 3 ? 0 : sedimentSelector.SelectedIndex;

            double maxFlux = 1e-6;
            double minFlux = -1e-6;
            foreach (var fluxes in sedimentFluxes)
            {
                for (int i = 0; i < gridPoints; i++)
                {
                    double flux = fluxes[selectedSediment * gridPoints + i];
                    maxFlux = Math.Max(maxFlux, flux);
                    minFlux = Math.Min(minFlux, flux);
                }
            }
            double fluxRange = Math.Max(maxFlux, -minFlux);

            for (int t = 0; t < sedimentFluxes.Count; t++)
            {
                float x = t * panelWidth / (float)maxTimeSteps;
                float width = panelWidth / (float)maxTimeSteps;
                for (int i = 0; i < gridPoints; i++)
                {
                    float y = i * panelHeight / (float)gridPoints;
                    float height = panelHeight / (float)gridPoints;
                    double flux = sedimentFluxes[t][selectedSediment * gridPoints + i];
                    int intensity = fluxRange > 0 ? (int)(255 * flux / fluxRange) : 0;
                    intensity = Math.Min(Math.Max(intensity, -255), 255);
                    Color color = intensity >= 0
                        ? Color.FromArgb(intensity, 0, 0, 255 - intensity)
                        : Color.FromArgb(-intensity, 255 - (-intensity), 0, 0);
                    using (Brush brush = new SolidBrush(color))
                    {
                        g.FillRectangle(brush, x, y, width, height);
                    }
                }
            }

            g.DrawLine(Pens.Black, 0, panelHeight - 10, panelWidth, panelHeight - 10);
            g.DrawLine(Pens.Black, 10, 0, 10, panelHeight);
            g.DrawString("Time (s)", new Font("Verdana", 10), Brushes.Black, panelWidth - 120, panelHeight - 30);
            g.DrawString("Distance (m)", new Font("Verdana", 10), Brushes.Black, 10, 10);
            g.DrawString($"{sedimentTypes[selectedSediment]} Flux", new Font("Verdana", 10), Brushes.Black, 10, 30);

            for (int i = 0; i < 100; i++)
            {
                int intensity = (i - 50) * 255 / 50;
                Color color = intensity >= 0
                    ? Color.FromArgb(intensity, 0, 0, 255 - intensity)
                    : Color.FromArgb(-intensity, 255 - (-intensity), 0, 0);
                using (Brush brush = new SolidBrush(color))
                {
                    g.FillRectangle(brush, panelWidth - 30, panelHeight - 110 + i, 20, 1);
                }
            }
            g.DrawString($"Flux (kg/m²s)\n[{minFlux:F6}, {maxFlux:F6}]", new Font("Verdana", 10), Brushes.Black, panelWidth - 70, panelHeight - 120);
        }

        private void BudgetMapPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            g.Clear(Color.White);

            int panelWidth = budgetMapPanel.Width;
            int panelHeight = budgetMapPanel.Height;
            int selectedSediment = sedimentSelector.SelectedIndex;
            bool showAllSediments = selectedSediment == 3;
            int numRows = showAllSediments ? sedimentTypes.Length : 1;
            float rowHeight = panelHeight / (float)numRows;

            if (sedimentFluxes.Count == 0) return;
            double[] latestFluxes = sedimentFluxes[sedimentFluxes.Count - 1];

            double maxFlux = 1e-6;
            double minFlux = -1e-6;
            double maxVelocity = 0.01;
            double maxTransport = 0.01;
            for (int f = 0; f < sedimentTypes.Length; f++)
            {
                if (!showAllSediments && f != selectedSediment) continue;
                for (int i = 0; i < gridPoints; i++)
                {
                    double flux = latestFluxes[f * gridPoints + i];
                    maxFlux = Math.Max(maxFlux, flux);
                    minFlux = Math.Min(minFlux, flux);
                    maxVelocity = Math.Max(maxVelocity, Math.Abs(velocityField[i]));
                    maxTransport = Math.Max(maxTransport, Math.Abs(velocityField[i] * sedimentConcentrations[f][i]));
                }
            }
            double fluxRange = Math.Max(maxFlux, -minFlux);

            for (int f = 0; f < sedimentTypes.Length; f++)
            {
                if (!showAllSediments && f != selectedSediment) continue;
                int row = showAllSediments ? f : 0;
                float yOffset = row * rowHeight;
                for (int i = 0; i < gridPoints; i++)
                {
                    float x = i * panelWidth / (float)gridPoints;
                    float width = panelWidth / (float)gridPoints;
                    double flux = latestFluxes[f * gridPoints + i];
                    Color color;
                    if (Math.Abs(flux) < fluxRange * 0.1)
                        color = Color.FromArgb(0, 255, 0);
                    else if (flux > 0)
                    {
                        int intensity = fluxRange > 0 ? (int)(255 * flux / fluxRange) : 0;
                        color = Color.FromArgb(intensity, 0, 0, 255 - intensity);
                    }
                    else
                    {
                        int intensity = fluxRange > 0 ? (int)(255 * -flux / fluxRange) : 0;
                        color = Color.FromArgb(intensity, 255 - intensity, 0, 0);
                    }
                    using (Brush brush = new SolidBrush(color))
                    {
                        g.FillRectangle(brush, x, yOffset, width, rowHeight);
                    }
                }
            }

            using (Pen velocityPen = new Pen(Color.Black, 2))
            using (Pen transportPen = new Pen(Color.Purple, 2))
            {
                for (int f = 0; f < sedimentTypes.Length; f++)
                {
                    if (!showAllSediments && f != selectedSediment) continue;
                    int row = showAllSediments ? f : 0;
                    float yOffset = row * rowHeight + rowHeight / 2;
                    for (int i = 0; i < gridPoints; i += 10)
                    {
                        float x = i * panelWidth / (float)gridPoints;
                        float vectorLength = 30.0f;
                        float velScale = (float)(vectorLength * velocityField[i] / maxVelocity);
                        g.DrawLine(velocityPen, x, yOffset, x + velScale, yOffset);
                        g.DrawLine(velocityPen, x + velScale, yOffset - 5, x + velScale, yOffset + 5);
                        double transport = velocityField[i] * sedimentConcentrations[f][i];
                        float transScale = (float)(vectorLength * transport / maxTransport);
                        g.DrawLine(transportPen, x, yOffset + 10, x + transScale, yOffset + 10);
                        g.DrawLine(transportPen, x + transScale, yOffset + 5, x + transScale, yOffset + 15);
                    }
                }
            }

            g.DrawLine(Pens.Black, 0, panelHeight - 10, panelWidth, panelHeight - 10);
            g.DrawLine(Pens.Black, 10, 0, 10, panelHeight);
            g.DrawString("Distance (m)", new Font("Verdana", 10), Brushes.Black, panelWidth - 120, panelHeight - 30);
            g.DrawString("Sediment Type", new Font("Verdana", 10), Brushes.Black, 10, 10);
            if (showAllSediments)
            {
                for (int f = 0; f < sedimentTypes.Length; f++)
                {
                    g.DrawString(sedimentTypes[f], new Font("Verdana", 10), Brushes.Black, 10, f * rowHeight + 5);
                }
            }
            else
            {
                g.DrawString(sedimentTypes[selectedSediment], new Font("Verdana", 10), Brushes.Black, 10, 5);
            }

            g.FillRectangle(Brushes.Blue, panelWidth - 70, 20, 15, 15);
            g.DrawString("Deposition", new Font("Verdana", 10), Brushes.Black, panelWidth - 50, 20);
            g.FillRectangle(Brushes.Red, panelWidth - 70, 40, 15, 15);
            g.DrawString("Erosion", new Font("Verdana", 10), Brushes.Black, panelWidth - 50, 40);
            g.FillRectangle(Brushes.Green, panelWidth - 70, 60, 15, 15);
            g.DrawString("Transport", new Font("Verdana", 10), Brushes.Black, panelWidth - 50, 60);
            g.DrawLine(new Pen(Color.Black, 2), panelWidth - 70, 80, panelWidth - 50, 80);
            g.DrawString("Velocity", new Font("Verdana", 10), Brushes.Black, panelWidth - 50, 80);
            g.DrawLine(new Pen(Color.Purple, 2), panelWidth - 70, 100, panelWidth - 50, 100);
            g.DrawString("Transport", new Font("Verdana", 10), Brushes.Black, panelWidth - 50, 100);

            for (int i = 0; i < 100; i++)
            {
                int intensity = (i - 50) * 255 / 50;
                Color color = intensity >= 0
                    ? Color.FromArgb(intensity, 0, 0, 255 - intensity)
                    : Color.FromArgb(-intensity, 255 - (-intensity), 0, 0);
                using (Brush brush = new SolidBrush(color))
                {
                    g.FillRectangle(brush, panelWidth - 30, panelHeight - 110 + i, 20, 1);
                }
            }
            g.DrawString($"Flux (kg/m²s)\n[{minFlux:F6}, {maxFlux:F6}]", new Font("Verdana", 10), Brushes.Black, panelWidth - 70, panelHeight - 120);
        }

        private void UpdateButtonStates()
        {
            startSimButton.Enabled = !isRunning;
            pauseSimButton.Enabled = isRunning;
            resetSimButton.Enabled = true;
        }

        public static void ShowMFSTWindow()
        {
            using (var window = new MultiFracSedimentTransp())
            {
                window.ShowDialog();
            }
        }
    }
}