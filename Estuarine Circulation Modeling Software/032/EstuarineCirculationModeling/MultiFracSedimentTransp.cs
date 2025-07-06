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
        private TabPage correlationTab;
        private Panel visualizationPanel;
        private Panel hovmoellerPanel;
        private Panel budgetMapPanel;
        private Panel correlationPanel;
        private TextBox sandConcentrationTextBox;
        private TextBox siltConcentrationTextBox;
        private TextBox clayConcentrationTextBox;
        private TextBox tidalAmplitudeTextBox;
        private TextBox tidalPeriodTextBox;
        private TextBox riverInflowTextBox;
        private TextBox salinityDiffusionTextBox;
        private ComboBox sedimentSelector;
        private CheckBox showHovmoellerCheckBox;
        private CheckBox showBudgetMapCheckBox;
        private CheckBox showVectorsCheckBox;
        private CheckBox showCorrelationCheckBox;
        private Label settlingVelocitiesLabel;
        private Button startSimButton;
        private Button pauseSimButton;
        private Button resetSimButton;
        private TextBox outputTextBox;
        private Timer simulationTimer;
        private bool isRunning = false;
        private int lastSelectedTabIndex = 0; // Persist selected tab
        private double[][] sedimentConcentrations; // [fraction][gridPoint]
        private double[] salinityField; // ppt
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
        private double tidalPeriod = 43200.0; // s (12 hours)
        private double riverInflow = 10.0; // m³/s
        private double salinityDiffusion = 0.1; // m²/s
        private const int maxTimeSteps = 432000; // 43200 s / 0.1 s

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
                Text = "43200"
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
                Text = "10.0"
            };
            riverInflowTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(riverInflowTextBox, "River inflow rate (0-100 m³/s)");
            var salinityDiffusionLabel = new Label
            {
                Location = new Point(10, 380),
                Size = new Size(280, 20),
                Text = "Salinity Diffusion (m²/s):"
            };
            salinityDiffusionTextBox = new TextBox
            {
                Location = new Point(10, 400),
                Size = new Size(260, 25),
                Text = "0.1"
            };
            salinityDiffusionTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(salinityDiffusionTextBox, "Salinity diffusion coefficient (0.01-1.0 m²/s)");

            // Visualization Options Group
            var vizGroupLabel = new Label
            {
                Location = new Point(10, 440),
                Size = new Size(280, 20),
                Text = "Visualization Options",
                Font = new Font("Verdana", 10F, FontStyle.Bold)
            };
            var sedimentSelectorLabel = new Label
            {
                Location = new Point(10, 470),
                Size = new Size(280, 20),
                Text = "Sediment for Maps:"
            };
            sedimentSelector = new ComboBox
            {
                Location = new Point(10, 490),
                Size = new Size(260, 25),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            sedimentSelector.Items.AddRange(new object[] { "Sand", "Silt", "Clay", "All" });
            sedimentSelector.SelectedIndex = 0;
            sedimentSelector.SelectedIndexChanged += (s, e) => { hovmoellerPanel.Invalidate(); budgetMapPanel.Invalidate(); correlationPanel.Invalidate(); };
            sedimentSelector.MouseHover += (s, e) => new ToolTip().SetToolTip(sedimentSelector, "Select sediment type for Hovmöller, budget, and correlation diagrams");

            showHovmoellerCheckBox = new CheckBox
            {
                Location = new Point(10, 520),
                Size = new Size(280, 25),
                Text = "Show Hovmöller Diagram",
                Checked = true
            };
            showHovmoellerCheckBox.MouseHover += (s, e) => new ToolTip().SetToolTip(showHovmoellerCheckBox, "Toggle temporal sediment flux diagram");
            showBudgetMapCheckBox = new CheckBox
            {
                Location = new Point(10, 550),
                Size = new Size(280, 25),
                Text = "Show Sediment Budget Map",
                Checked = true
            };
            showBudgetMapCheckBox.MouseHover += (s, e) => new ToolTip().SetToolTip(showBudgetMapCheckBox, "Toggle spatial sediment budget map");
            showVectorsCheckBox = new CheckBox
            {
                Location = new Point(10, 580),
                Size = new Size(280, 25),
                Text = "Show Vectors in Budget Map",
                Checked = true
            };
            showVectorsCheckBox.MouseHover += (s, e) => new ToolTip().SetToolTip(showVectorsCheckBox, "Toggle velocity and transport vectors in budget map");
            showVectorsCheckBox.CheckedChanged += (s, e) => budgetMapPanel.Invalidate();
            showCorrelationCheckBox = new CheckBox
            {
                Location = new Point(10, 610),
                Size = new Size(280, 25),
                Text = "Show Correlation Diagram",
                Checked = true
            };
            showCorrelationCheckBox.MouseHover += (s, e) => new ToolTip().SetToolTip(showCorrelationCheckBox, "Toggle salinity-sediment correlation diagram");
            showCorrelationCheckBox.CheckedChanged += (s, e) => {
                correlationTab.Enabled = showCorrelationCheckBox.Checked;
                if (!showCorrelationCheckBox.Checked && tabControl.SelectedTab == correlationTab)
                    tabControl.SelectedTab = visualizationTab;
            };

            settlingVelocitiesLabel = new Label
            {
                Location = new Point(10, 640),
                Size = new Size(280, 80),
                Text = $"Settling Velocities (m/s):\nSand: {settlingVelocities[0]:F6}\nSilt: {settlingVelocities[1]:F6}\nClay: {settlingVelocities[2]:F6}"
            };

            // Control Buttons
            startSimButton = new Button
            {
                Location = new Point(10, 730),
                Size = new Size(120, 30),
                Text = "Start/Resume",
                FlatStyle = FlatStyle.Flat
            };
            startSimButton.Click += StartSimButton_Click;
            startSimButton.MouseHover += (s, e) => new ToolTip().SetToolTip(startSimButton, "Start or resume the simulation");
            pauseSimButton = new Button
            {
                Location = new Point(140, 730),
                Size = new Size(120, 30),
                Text = "Pause",
                FlatStyle = FlatStyle.Flat,
                Enabled = false
            };
            pauseSimButton.Click += PauseSimButton_Click;
            pauseSimButton.MouseHover += (s, e) => new ToolTip().SetToolTip(pauseSimButton, "Pause the simulation");
            resetSimButton = new Button
            {
                Location = new Point(10, 770),
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
                Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Right | AnchorStyles.Bottom,
                Appearance = TabAppearance.Normal,
                SelectedIndex = lastSelectedTabIndex
            };
            tabControl.SelectedIndexChanged += (s, e) => lastSelectedTabIndex = tabControl.SelectedIndex;

            visualizationTab = new TabPage { Text = "Profiles" };
            new ToolTip().SetToolTip(visualizationTab, "Concentration, velocity, and salinity profiles");
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
            new ToolTip().SetToolTip(hovmoellerTab, "Temporal-spatial sediment flux diagram");
            hovmoellerPanel = new Panel
            {
                Location = new Point(0, 0),
                Size = new Size(856, 666),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White,
                Visible = showHovmoellerCheckBox.Checked
            };
            hovmoellerPanel.Paint += HovmoellerPanel_Paint;
            showHovmoellerCheckBox.CheckedChanged += (s, e) => {
                hovmoellerTab.Enabled = showHovmoellerCheckBox.Checked;
                if (!showHovmoellerCheckBox.Checked && tabControl.SelectedTab == hovmoellerTab)
                    tabControl.SelectedTab = visualizationTab;
            };
            hovmoellerTab.Controls.Add(hovmoellerPanel);

            budgetMapTab = new TabPage { Text = "Budget Map" };
            new ToolTip().SetToolTip(budgetMapTab, "Spatial sediment budget with zones and vectors");
            budgetMapPanel = new Panel
            {
                Location = new Point(0, 0),
                Size = new Size(856, 666),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White,
                Visible = showBudgetMapCheckBox.Checked
            };
            budgetMapPanel.Paint += BudgetMapPanel_Paint;
            showBudgetMapCheckBox.CheckedChanged += (s, e) => {
                budgetMapTab.Enabled = showBudgetMapCheckBox.Checked;
                if (!showBudgetMapCheckBox.Checked && tabControl.SelectedTab == budgetMapTab)
                    tabControl.SelectedTab = visualizationTab;
            };
            budgetMapTab.Controls.Add(budgetMapPanel);

            correlationTab = new TabPage { Text = "Correlation" };
            new ToolTip().SetToolTip(correlationTab, "Salinity vs. sediment concentration correlation");
            correlationPanel = new Panel
            {
                Location = new Point(0, 0),
                Size = new Size(856, 666),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White,
                Visible = showCorrelationCheckBox.Checked
            };
            correlationPanel.Paint += CorrelationPanel_Paint;
            correlationTab.Controls.Add(correlationPanel);

            tabControl.TabPages.AddRange(new[] { visualizationTab, hovmoellerTab, budgetMapTab, correlationTab });

            controlPanel.Controls.AddRange(new Control[] {
                sedimentGroupLabel, sandLabel, sandConcentrationTextBox,
                siltLabel, siltConcentrationTextBox,
                clayLabel, clayConcentrationTextBox,
                flowGroupLabel, tidalAmplitudeLabel, tidalAmplitudeTextBox,
                tidalPeriodLabel, tidalPeriodTextBox,
                riverInflowLabel, riverInflowTextBox,
                salinityDiffusionLabel, salinityDiffusionTextBox,
                vizGroupLabel, sedimentSelectorLabel, sedimentSelector,
                showHovmoellerCheckBox, showBudgetMapCheckBox, showVectorsCheckBox, showCorrelationCheckBox,
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
            salinityField = new double[gridPoints];
            sedimentFluxes = new List<double[]>();
            for (int i = 0; i < gridPoints; i++)
            {
                positions[i] = i * estuaryLength / (gridPoints - 1);
                velocityField[i] = riverInflow / estuaryDepth;
                salinityField[i] = 35.0 * (1.0 - Math.Exp(-positions[i] / estuaryLength)); // Exponential initial profile
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
                salinityDiffusion = double.Parse(salinityDiffusionTextBox.Text.Trim());

                if (sand < 0 || silt < 0 || clay < 0)
                    throw new Exception("Concentrations must be non-negative.");
                if (tidalAmplitude < 0 || tidalAmplitude > 10.0)
                    throw new Exception("Tidal amplitude must be between 0 and 10 m.");
                if (tidalPeriod < 600.0 || tidalPeriod > 86400.0)
                    throw new Exception("Tidal period must be between 600 and 86400 s.");
                if (riverInflow < 0 || riverInflow > 100.0)
                    throw new Exception("River inflow must be between 0 and 100 m³/s.");
                if (salinityDiffusion < 0.01 || salinityDiffusion > 1.0)
                    throw new Exception("Salinity diffusion must be between 0.01 and 1.0 m²/s.");

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
                    correlationPanel.Invalidate();
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
            for (int i = 0; i < gridPoints; i++)
                salinityField[i] = 35.0 * (1.0 - Math.Exp(-positions[i] / estuaryLength));
            outputTextBox.Clear();
            outputTextBox.AppendText("Simulation reset.\r\n");
            visualizationPanel.Invalidate();
            hovmoellerPanel.Invalidate();
            budgetMapPanel.Invalidate();
            correlationPanel.Invalidate();
            tabControl.SelectedIndex = lastSelectedTabIndex;
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
                double distanceFactor = Math.Exp(-x / estuaryLength); // Exponential decay for tidal influence
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

            // Update sediment concentrations
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

            // Update salinity
            double[] newSalinity = new double[gridPoints];
            double tidalPhase = 2 * Math.PI * time / tidalPeriod;
            double intrusionFactor = 0.5 * (1 + Math.Cos(tidalPhase)); // Enhances salinity intrusion during flood
            for (int i = 1; i < gridPoints - 1; i++)
            {
                double advection = velocityField[i] >= 0
                    ? -velocityField[i] * (salinityField[i] - salinityField[i - 1]) / dx
                    : -velocityField[i] * (salinityField[i + 1] - salinityField[i]) / dx;
                double diffusion = salinityDiffusion * (salinityField[i + 1] - 2 * salinityField[i] + salinityField[i - 1]) / (dx * dx);
                newSalinity[i] = salinityField[i] + timeStep * (advection + diffusion);
                if (newSalinity[i] < 0) newSalinity[i] = 0;
                if (newSalinity[i] > 35) newSalinity[i] = 35;
            }
            newSalinity[0] = 0.0; // River boundary
            newSalinity[gridPoints - 1] = 35.0 * (0.8 + 0.2 * intrusionFactor); // Dynamic ocean boundary
            salinityField = newSalinity;

            sedimentFluxes.Add(currentFluxes);
            if (sedimentFluxes.Count > maxTimeSteps)
                sedimentFluxes.RemoveAt(0);

            time += timeStep;
            visualizationPanel.Invalidate();
            hovmoellerPanel.Invalidate();
            budgetMapPanel.Invalidate();
            correlationPanel.Invalidate();
            UpdateOutput();
        }

        private void UpdateOutput()
        {
            double tidalPhase = (2 * Math.PI * time / tidalPeriod) % (2 * Math.PI);
            string phaseStr = tidalPhase < Math.PI ? "Flood" : "Ebb";
            double maxVelocity = velocityField[0];
            double minVelocity = velocityField[gridPoints - 1];
            double maxSalinity = salinityField[0];
            double minSalinity = salinityField[0];
            double maxETMConc = 0.0;
            int etmIndex = -1;
            for (int i = 0; i < gridPoints; i++)
            {
                maxSalinity = Math.Max(maxSalinity, salinityField[i]);
                minSalinity = Math.Min(minSalinity, salinityField[i]);
                if (salinityField[i] >= 5.0 && salinityField[i] <= 15.0)
                {
                    for (int f = 0; f < sedimentTypes.Length; f++)
                    {
                        if (sedimentConcentrations[f][i] > maxETMConc)
                        {
                            maxETMConc = sedimentConcentrations[f][i];
                            etmIndex = i;
                        }
                    }
                }
            }
            string output = $"Time: {time:F2}s | Tidal Phase: {phaseStr} | Vel (mid): {velocityField[gridPoints / 2]:F6} m/s | Vel (range): [{minVelocity:F6}, {maxVelocity:F6}] m/s | Salinity (mid): {salinityField[gridPoints / 2]:F2} ppt | Salinity (range): [{minSalinity:F2}, {maxSalinity:F2}] ppt | ETM Conc: {maxETMConc:F6} kg/m³ @ x={positions[etmIndex]:F0}m\r\n";
            outputTextBox.AppendText(output);
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

            double minSalinity = salinityField[0];
            double maxSalinity = salinityField[0];
            for (int i = 0; i < gridPoints; i++)
            {
                minSalinity = Math.Min(minSalinity, salinityField[i]);
                maxSalinity = Math.Max(maxSalinity, salinityField[i]);
            }
            double salinityRange = maxSalinity - minSalinity;
            if (salinityRange < 1e-6) salinityRange = 1.0; // Prevent division by zero
            using (Pen salinityPen = new Pen(Color.FromArgb(128, 0, 128, 255), 2))
            {
                for (int i = 0; i < gridPoints - 1; i++)
                {
                    float x1 = i * panelWidth / (float)(gridPoints - 1);
                    float x2 = (i + 1) * panelWidth / (float)(gridPoints - 1);
                    float y1 = panelHeight - (float)((salinityField[i] - minSalinity) / salinityRange * panelHeight * 0.4 + panelHeight * 0.1);
                    float y2 = panelHeight - (float)((salinityField[i + 1] - minSalinity) / salinityRange * panelHeight * 0.4 + panelHeight * 0.1);
                    g.DrawLine(salinityPen, x1, y1, x2, y2);
                }
            }

            g.DrawLine(Pens.Black, 0, panelHeight - 10, panelWidth, panelHeight - 10);
            g.DrawLine(Pens.Black, 10, 0, 10, panelHeight);
            g.DrawString("Distance (m)", new Font("Verdana", 11), Brushes.Black, panelWidth - 120, panelHeight - 30);
            g.DrawString("Concentration (kg/m³)", new Font("Verdana", 11), Brushes.Black, 10, 10);
            g.DrawString("Velocity (m/s)", new Font("Verdana", 11), Brushes.Black, 10, panelHeight / 2 - 20);
            g.DrawString($"Salinity (ppt) [{minSalinity:F2}, {maxSalinity:F2}]", new Font("Verdana", 11), Brushes.Black, 10, panelHeight - 50);

            for (int f = 0; f < sedimentTypes.Length; f++)
            {
                g.FillRectangle(new SolidBrush(colors[f]), 20, 20 + f * 25, 15, 15);
                g.DrawString(sedimentTypes[f], new Font("Verdana", 11), Brushes.Black, 40, 20 + f * 25);
            }
            g.FillRectangle(new SolidBrush(Color.Black), 20, 95, 15, 15);
            g.DrawString("Velocity", new Font("Verdana", 11), Brushes.Black, 40, 95);
            g.FillRectangle(new SolidBrush(Color.FromArgb(128, 0, 128, 255)), 20, 120, 15, 15);
            g.DrawString("Salinity", new Font("Verdana", 11), Brushes.Black, 40, 120);

            double tidalPhase = (2 * Math.PI * time / tidalPeriod) % (2 * Math.PI);
            string phaseStr = tidalPhase < Math.PI ? "Flood" : "Ebb";
            g.DrawString($"Time: {time:F2} s | Phase: {phaseStr}", new Font("Verdana", 11), Brushes.Black, panelWidth - 180, 10);
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
            g.DrawString("Time (s)", new Font("Verdana", 11), Brushes.Black, panelWidth - 120, panelHeight - 30);
            g.DrawString("Distance (m)", new Font("Verdana", 11), Brushes.Black, 10, 10);
            g.DrawString($"{sedimentTypes[selectedSediment]} Flux", new Font("Verdana", 11), Brushes.Black, 10, 30);

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
            g.DrawString($"Flux (kg/m²s)\n[{minFlux:F6}, {maxFlux:F6}]", new Font("Verdana", 11), Brushes.Black, panelWidth - 70, panelHeight - 120);
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

            if (showVectorsCheckBox.Checked)
            {
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
            }

            g.DrawLine(Pens.Black, 0, panelHeight - 10, panelWidth, panelHeight - 10);
            g.DrawLine(Pens.Black, 10, 0, 10, panelHeight);
            g.DrawString("Distance (m)", new Font("Verdana", 11), Brushes.Black, panelWidth - 120, panelHeight - 30);
            g.DrawString("Sediment Type", new Font("Verdana", 11), Brushes.Black, 10, 10);
            if (showAllSediments)
            {
                for (int f = 0; f < sedimentTypes.Length; f++)
                {
                    g.DrawString(sedimentTypes[f], new Font("Verdana", 11), Brushes.Black, 10, f * rowHeight + 5);
                }
            }
            else
            {
                g.DrawString(sedimentTypes[selectedSediment], new Font("Verdana", 11), Brushes.Black, 10, 5);
            }

            g.FillRectangle(Brushes.Blue, panelWidth - 70, 20, 15, 15);
            g.DrawString("Deposition", new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 20);
            g.FillRectangle(Brushes.Red, panelWidth - 70, 40, 15, 15);
            g.DrawString("Erosion", new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 40);
            g.FillRectangle(Brushes.Green, panelWidth - 70, 60, 15, 15);
            g.DrawString("Transport", new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 60);
            if (showVectorsCheckBox.Checked)
            {
                g.DrawLine(new Pen(Color.Black, 2), panelWidth - 70, 80, panelWidth - 50, 80);
                g.DrawString("Velocity", new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 80);
                g.DrawLine(new Pen(Color.Purple, 2), panelWidth - 70, 100, panelWidth - 50, 100);
                g.DrawString("Transport", new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 100);
            }

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
            g.DrawString($"Flux (kg/m²s)\n[{minFlux:F6}, {maxFlux:F6}]", new Font("Verdana", 11), Brushes.Black, panelWidth - 70, panelHeight - 120);
        }

        private void CorrelationPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            g.Clear(Color.White);

            int panelWidth = correlationPanel.Width;
            int panelHeight = correlationPanel.Height;
            int selectedSediment = sedimentSelector.SelectedIndex;
            bool showAllSediments = selectedSediment == 3;
            double maxConcentration = 0.2; // kg/m³
            double maxSalinity = 35.0; // ppt
            double minSalinity = 0.0;
            for (int i = 0; i < gridPoints; i++)
            {
                minSalinity = Math.Min(minSalinity, salinityField[i]);
                maxSalinity = Math.Max(maxSalinity, salinityField[i]);
            }
            double salinityRange = maxSalinity - minSalinity;
            if (salinityRange < 1e-6) salinityRange = 35.0; // Default range

            if (showAllSediments)
            {
                // Correlation matrix for all sediments
                int gridSize = panelWidth / 3;
                double maxETMConc = 0.0;
                for (int i = 0; i < gridPoints; i++)
                {
                    if (salinityField[i] >= 5.0 && salinityField[i] <= 15.0)
                    {
                        for (int f = 0; f < sedimentTypes.Length; f++)
                            maxETMConc = Math.Max(maxETMConc, sedimentConcentrations[f][i]);
                    }
                }
                for (int f = 0; f < sedimentTypes.Length; f++)
                {
                    // Scatter plots: sediment vs. salinity
                    int row = f;
                    float plotX = gridSize;
                    float plotY = row * gridSize;
                    float plotWidth = gridSize - 20;
                    float plotHeight = gridSize - 20;
                    g.DrawRectangle(Pens.Black, plotX, plotY, plotWidth, plotHeight);

                    // ETM salinity band (5-15 ppt)
                    float bandX1 = plotX + (float)((5.0 - minSalinity) / salinityRange * plotWidth);
                    float bandX2 = plotX + (float)((15.0 - minSalinity) / salinityRange * plotWidth);
                    g.FillRectangle(new SolidBrush(Color.FromArgb(50, 0, 255, 0)), bandX1, plotY, bandX2 - bandX1, plotHeight);

                    for (int i = 0; i < gridPoints; i++)
                    {
                        float x = plotX + (float)((salinityField[i] - minSalinity) / salinityRange * plotWidth);
                        float y = plotY + plotHeight - (float)(sedimentConcentrations[f][i] / maxConcentration * plotHeight);
                        int size = (int)(5 + 10 * sedimentConcentrations[f][i] / maxConcentration);
                        Color color = f == 0 ? Color.Red : f == 1 ? Color.Blue : Color.Green;
                        int alpha = salinityField[i] >= 5.0 && salinityField[i] <= 15.0 ? 255 : 128;
                        using (Brush brush = new SolidBrush(Color.FromArgb(alpha, color)))
                        {
                            g.FillEllipse(brush, x - size / 2, y - size / 2, size, size);
                        }
                    }

                    g.DrawString($"{sedimentTypes[f]} vs. Salinity", new Font("Verdana", 10), Brushes.Black, plotX + 5, plotY + 5);
                    g.DrawString("Salinity (ppt)", new Font("Verdana", 10), Brushes.Black, plotX, plotY + plotHeight + 5);
                    g.DrawString("Conc (kg/m³)", new Font("Verdana", 10), Brushes.Black, plotX - 50, plotY + plotHeight / 2);
                }

                // ETM concentration line
                float etmY = panelHeight - (float)(maxETMConc / maxConcentration * (panelHeight - 100)) + 50;
                g.DrawLine(new Pen(Color.Black, 1) { DashStyle = System.Drawing.Drawing2D.DashStyle.Dash }, 50, etmY, panelWidth - 50, etmY);
                g.DrawString($"Max ETM Conc: {maxETMConc:F6} kg/m³", new Font("Verdana", 10), Brushes.Black, panelWidth - 200, etmY - 20);

                // Labels and legend
                g.DrawString("Salinity-Sediment Correlation Matrix", new Font("Verdana", 11), Brushes.Black, 10, 10);
                g.FillRectangle(new SolidBrush(Color.FromArgb(50, 0, 255, 0)), panelWidth - 70, 20, 15, 15);
                g.DrawString("ETM Zone (5-15 ppt)", new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 20);
                for (int f = 0; f < sedimentTypes.Length; f++)
                {
                    Color color = f == 0 ? Color.Red : f == 1 ? Color.Blue : Color.Green;
                    g.FillEllipse(new SolidBrush(color), panelWidth - 70, 40 + f * 20, 10, 10);
                    g.DrawString(sedimentTypes[f], new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 40 + f * 20);
                }
            }
            else
            {
                // Single sediment scatter plot
                float plotX = 50;
                float plotY = 50;
                float plotWidth = panelWidth - 100;
                float plotHeight = panelHeight - 100;
                g.DrawRectangle(Pens.Black, plotX, plotY, plotWidth, plotHeight);

                // ETM salinity band (5-15 ppt)
                float bandX1 = plotX + (float)((5.0 - minSalinity) / salinityRange * plotWidth);
                float bandX2 = plotX + (float)((15.0 - minSalinity) / salinityRange * plotWidth);
                g.FillRectangle(new SolidBrush(Color.FromArgb(50, 0, 255, 0)), bandX1, plotY, bandX2 - bandX1, plotHeight);

                Color color = selectedSediment == 0 ? Color.Red : selectedSediment == 1 ? Color.Blue : Color.Green;
                double maxETMConc = 0.0;
                for (int i = 0; i < gridPoints; i++)
                {
                    if (salinityField[i] >= 5.0 && salinityField[i] <= 15.0)
                        maxETMConc = Math.Max(maxETMConc, sedimentConcentrations[selectedSediment][i]);
                    float x = plotX + (float)((salinityField[i] - minSalinity) / salinityRange * plotWidth);
                    float y = plotY + plotHeight - (float)(sedimentConcentrations[selectedSediment][i] / maxConcentration * plotHeight);
                    int size = (int)(5 + 10 * sedimentConcentrations[selectedSediment][i] / maxConcentration);
                    int alpha = salinityField[i] >= 5.0 && salinityField[i] <= 15.0 ? 255 : 128;
                    using (Brush brush = new SolidBrush(Color.FromArgb(alpha, color)))
                    {
                        g.FillEllipse(brush, x - size / 2, y - size / 2, size, size);
                    }
                }

                // ETM concentration line
                float etmY = plotY + plotHeight - (float)(maxETMConc / maxConcentration * plotHeight);
                g.DrawLine(new Pen(Color.Black, 1) { DashStyle = System.Drawing.Drawing2D.DashStyle.Dash }, plotX, etmY, plotX + plotWidth, etmY);
                g.DrawString($"Max ETM Conc: {maxETMConc:F6} kg/m³", new Font("Verdana", 10), Brushes.Black, plotX + plotWidth - 200, etmY - 20);

                g.DrawString($"{sedimentTypes[selectedSediment]} vs. Salinity", new Font("Verdana", 11), Brushes.Black, plotX + 5, plotY - 20);
                g.DrawString($"Salinity (ppt) [{minSalinity:F2}, {maxSalinity:F2}]", new Font("Verdana", 11), Brushes.Black, plotX, plotY + plotHeight + 5);
                g.DrawString("Concentration (kg/m³)", new Font("Verdana", 11), Brushes.Black, plotX - 50, plotY + plotHeight / 2);

                // Legend
                g.FillRectangle(new SolidBrush(Color.FromArgb(50, 0, 255, 0)), panelWidth - 70, 20, 15, 15);
                g.DrawString("ETM Zone (5-15 ppt)", new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 20);
                g.FillEllipse(new SolidBrush(color), panelWidth - 70, 40, 10, 10);
                g.DrawString(sedimentTypes[selectedSediment], new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 40);
            }
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