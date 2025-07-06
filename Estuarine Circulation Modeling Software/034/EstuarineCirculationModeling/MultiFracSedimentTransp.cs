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
        private TextBox flocculationCoefficientTextBox;
        private TextBox minFlocFactorTextBox;
        private ComboBox sedimentSelector;
        private CheckBox showFluxCheckBox;
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
        private int lastSelectedTabIndex = 0;
        private double[][] sedimentConcentrations; // [fraction][gridPoint]
        private double[] salinityField; // ppt
        private double[] shearRate; // s^-1
        private double[] flocSizes; // m (for clay)
        private readonly string[] sedimentTypes = { "Sand", "Silt", "Clay" };
        private readonly double[] grainSizes = { 0.0002, 0.00002, 0.000002 }; // m
        private double[] settlingVelocities; // m/s
        private readonly double[] diffusionCoefficients = { 0.001, 0.0005, 0.0002 }; // m²/s
        private readonly double[] criticalShearStresses = { 0.2, 0.1, 0.05 }; // N/m²
        private double[] positions;
        private double[] velocityField; // m/s
        private List<double[]> sedimentFluxes; // [time][fraction * gridPoint]
        private double time = 0.0;
        private double estuaryLength = 10000.0; // m
        private double estuaryDepth = 10.0; // m
        private double timeStep = 0.1; // s
        private int gridPoints = 100;
        private double tidalAmplitude = 1.0; // m
        private double tidalPeriod = 43200.0; // s (12 hours)
        private double riverInflow = 10.0; // m³/s
        private double salinityDiffusion = 0.1; // m²/s
        private double flocculationCoefficient = 0.5; // m³/kg
        private double minFlocFactor = 1.0; // Minimum flocculation factor
        private const double flocSizeBase = 0.000002; // m (2 µm)
        private const double optimalSalinity = 10.0; // ppt
        private const double salinityScale = 5.0; // ppt
        private const double optimalShearRate = 10.0; // s^-1
        private const double breakupShearRate = 200.0; // s^-1
        private const double referenceConcentration = 0.1; // kg/m³
        private const double concentrationExponent = 0.5;
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
            double waterViscosity = 1.0e-6; // m²/s
            double sedimentDensity = 2650.0; // kg/m³
            double waterDensity = 1000.0; // kg/m³
            double reductionFactor = 0.001;

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
                Anchor = AnchorStyles.Top | AnchorStyles.Left | AnchorStyles.Bottom,
                AutoScroll = true
            };

            var innerControlPanel = new Panel
            {
                Location = new Point(0, 0),
                Size = new Size(280, 900), // Increased height to accommodate all controls
                AutoSize = false
            };

            // Sediment Inputs Group
            var sedimentGroupLabel = new Label
            {
                Location = new Point(10, 10),
                Size = new Size(260, 20),
                Text = "Sediment Concentrations",
                Font = new Font("Verdana", 10F, FontStyle.Bold)
            };
            var sandLabel = new Label
            {
                Location = new Point(10, 40),
                Size = new Size(260, 20),
                Text = "Sand (kg/m³):"
            };
            sandConcentrationTextBox = new TextBox
            {
                Location = new Point(10, 60),
                Size = new Size(240, 25),
                Text = "0.1"
            };
            sandConcentrationTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(sandConcentrationTextBox, "Initial sand concentration (kg/m³)");
            var siltLabel = new Label
            {
                Location = new Point(10, 90),
                Size = new Size(260, 20),
                Text = "Silt (kg/m³):"
            };
            siltConcentrationTextBox = new TextBox
            {
                Location = new Point(10, 110),
                Size = new Size(240, 25),
                Text = "0.05"
            };
            siltConcentrationTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(siltConcentrationTextBox, "Initial silt concentration (kg/m³)");
            var clayLabel = new Label
            {
                Location = new Point(10, 140),
                Size = new Size(260, 20),
                Text = "Clay (kg/m³):"
            };
            clayConcentrationTextBox = new TextBox
            {
                Location = new Point(10, 160),
                Size = new Size(240, 25),
                Text = "0.02"
            };
            clayConcentrationTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(clayConcentrationTextBox, "Initial clay concentration (kg/m³)");

            // Flow Inputs Group
            var flowGroupLabel = new Label
            {
                Location = new Point(10, 200),
                Size = new Size(260, 20),
                Text = "Flow Parameters",
                Font = new Font("Verdana", 10F, FontStyle.Bold)
            };
            var tidalAmplitudeLabel = new Label
            {
                Location = new Point(10, 230),
                Size = new Size(260, 20),
                Text = "Tidal Amplitude (m):"
            };
            tidalAmplitudeTextBox = new TextBox
            {
                Location = new Point(10, 250),
                Size = new Size(240, 25),
                Text = "1.0"
            };
            tidalAmplitudeTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(tidalAmplitudeTextBox, "Tidal wave amplitude (0-10 m)");
            var tidalPeriodLabel = new Label
            {
                Location = new Point(10, 280),
                Size = new Size(260, 20),
                Text = "Tidal Period (s):"
            };
            tidalPeriodTextBox = new TextBox
            {
                Location = new Point(10, 300),
                Size = new Size(240, 25),
                Text = "43200"
            };
            tidalPeriodTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(tidalPeriodTextBox, "Tidal period (600-86400 s)");
            var riverInflowLabel = new Label
            {
                Location = new Point(10, 330),
                Size = new Size(260, 20),
                Text = "River Inflow (m³/s):"
            };
            riverInflowTextBox = new TextBox
            {
                Location = new Point(10, 350),
                Size = new Size(240, 25),
                Text = "10.0"
            };
            riverInflowTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(riverInflowTextBox, "River inflow rate (0-100 m³/s)");
            var salinityDiffusionLabel = new Label
            {
                Location = new Point(10, 380),
                Size = new Size(260, 20),
                Text = "Salinity Diffusion (m²/s):"
            };
            salinityDiffusionTextBox = new TextBox
            {
                Location = new Point(10, 400),
                Size = new Size(240, 25),
                Text = "0.1"
            };
            salinityDiffusionTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(salinityDiffusionTextBox, "Salinity diffusion coefficient (0.01-1.0 m²/s)");
            var flocculationCoefficientLabel = new Label
            {
                Location = new Point(10, 430),
                Size = new Size(260, 20),
                Text = "Flocculation Coeff (m³/kg):"
            };
            flocculationCoefficientTextBox = new TextBox
            {
                Location = new Point(10, 450),
                Size = new Size(240, 25),
                Text = "0.8"
            };
            flocculationCoefficientTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(flocculationCoefficientTextBox, "Flocculation aggregation coefficient (0.01-1.0 m³/kg)");
            var minFlocFactorLabel = new Label
            {
                Location = new Point(10, 480),
                Size = new Size(260, 20),
                Text = "Min Floc Factor:"
            };
            minFlocFactorTextBox = new TextBox
            {
                Location = new Point(10, 500),
                Size = new Size(240, 25),
                Text = "2.0"
            };
            minFlocFactorTextBox.MouseHover += (s, e) => new ToolTip().SetToolTip(minFlocFactorTextBox, "Minimum flocculation factor (0.1-10.0)");

            // Visualization Options Group
            var vizGroupLabel = new Label
            {
                Location = new Point(10, 540),
                Size = new Size(260, 20),
                Text = "Visualization Options",
                Font = new Font("Verdana", 10F, FontStyle.Bold)
            };
            var sedimentSelectorLabel = new Label
            {
                Location = new Point(10, 570),
                Size = new Size(260, 20),
                Text = "Sediment for Maps:"
            };
            sedimentSelector = new ComboBox
            {
                Location = new Point(10, 590),
                Size = new Size(240, 25),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            sedimentSelector.Items.AddRange(new object[] { "Sand", "Silt", "Clay", "All" });
            sedimentSelector.SelectedIndex = 0;
            sedimentSelector.SelectedIndexChanged += (s, e) => { hovmoellerPanel.Invalidate(); budgetMapPanel.Invalidate(); correlationPanel.Invalidate(); };
            sedimentSelector.MouseHover += (s, e) => new ToolTip().SetToolTip(sedimentSelector, "Select sediment type for Hovmöller, budget, and correlation diagrams");

            showFluxCheckBox = new CheckBox
            {
                Location = new Point(10, 620),
                Size = new Size(260, 25),
                Text = "Enable Flux Curve",
                Checked = true
            };
            showFluxCheckBox.MouseHover += (s, e) => new ToolTip().SetToolTip(showFluxCheckBox, "Toggle temporal sediment flux diagram");
            showBudgetMapCheckBox = new CheckBox
            {
                Location = new Point(10, 650),
                Size = new Size(260, 25),
                Text = "Show Sediment Budget Map",
                Checked = true
            };
            showBudgetMapCheckBox.MouseHover += (s, e) => new ToolTip().SetToolTip(showBudgetMapCheckBox, "Toggle spatial sediment budget map");
            showVectorsCheckBox = new CheckBox
            {
                Location = new Point(10, 680),
                Size = new Size(260, 25),
                Text = "Show Vectors in Budget Map",
                Checked = true
            };
            showVectorsCheckBox.MouseHover += (s, e) => new ToolTip().SetToolTip(showVectorsCheckBox, "Toggle velocity and transport vectors in budget map");
            showVectorsCheckBox.CheckedChanged += (s, e) => budgetMapPanel.Invalidate();
            showCorrelationCheckBox = new CheckBox
            {
                Location = new Point(10, 710),
                Size = new Size(260, 25),
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
                Location = new Point(10, 740),
                Size = new Size(260, 80),
                Text = $"Base Settling Velocities (m/s):\nSand: {settlingVelocities[0]:F6}\nSilt: {settlingVelocities[1]:F6}\nClay: {settlingVelocities[2]:F6}"
            };

            // Control Buttons
            startSimButton = new Button
            {
                Location = new Point(10, 820),
                Size = new Size(110, 30),
                Text = "Start",
                FlatStyle = FlatStyle.Flat
            };
            startSimButton.Click += StartSimButton_Click;
            startSimButton.MouseHover += (s, e) => new ToolTip().SetToolTip(startSimButton, "Start or resume the simulation");
            pauseSimButton = new Button
            {
                Location = new Point(130, 820),
                Size = new Size(110, 30),
                Text = "Pause",
                FlatStyle = FlatStyle.Flat,
                Enabled = false
            };
            pauseSimButton.Click += PauseSimButton_Click;
            pauseSimButton.MouseHover += (s, e) => new ToolTip().SetToolTip(pauseSimButton, "Pause the simulation");
            resetSimButton = new Button
            {
                Location = new Point(10, 860),
                Size = new Size(110, 30),
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

            hovmoellerTab = new TabPage { Text = "Flux" };
            new ToolTip().SetToolTip(hovmoellerTab, "Temporal-spatial sediment flux diagram");
            hovmoellerPanel = new Panel
            {
                Location = new Point(0, 0),
                Size = new Size(856, 666),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White,
                Visible = showFluxCheckBox.Checked
            };
            hovmoellerPanel.Paint += HovmoellerPanel_Paint;
            showFluxCheckBox.CheckedChanged += (s, e) => {
                hovmoellerTab.Enabled = showFluxCheckBox.Checked;
                if (!showFluxCheckBox.Checked && tabControl.SelectedTab == hovmoellerTab)
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
            new ToolTip().SetToolTip(correlationTab, "Salinity vs. sediment concentration and floc size");
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

            innerControlPanel.Controls.AddRange(new Control[] {
                sedimentGroupLabel, sandLabel, sandConcentrationTextBox,
                siltLabel, siltConcentrationTextBox,
                clayLabel, clayConcentrationTextBox,
                flowGroupLabel, tidalAmplitudeLabel, tidalAmplitudeTextBox,
                tidalPeriodLabel, tidalPeriodTextBox,
                riverInflowLabel, riverInflowTextBox,
                salinityDiffusionLabel, salinityDiffusionTextBox,
                flocculationCoefficientLabel, flocculationCoefficientTextBox,
                minFlocFactorLabel, minFlocFactorTextBox,
                vizGroupLabel, sedimentSelectorLabel, sedimentSelector,
                showFluxCheckBox, showBudgetMapCheckBox, showVectorsCheckBox, showCorrelationCheckBox,
                settlingVelocitiesLabel,
                startSimButton, pauseSimButton, resetSimButton
            });

            controlPanel.Controls.Add(innerControlPanel);
            this.Controls.AddRange(new Control[] { controlPanel, outputTextBox, tabControl });
        }

        // The rest of the class remains unchanged
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
                flocculationCoefficient = double.Parse(flocculationCoefficientTextBox.Text.Trim());
                minFlocFactor = double.Parse(minFlocFactorTextBox.Text.Trim());

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
                if (flocculationCoefficient < 0.01 || flocculationCoefficient > 1.0)
                    throw new Exception("Flocculation coefficient must be between 0.01 and 1.0 m³/kg.");
                if (minFlocFactor < 0.1 || minFlocFactor > 10.0)
                    throw new Exception("Minimum flocculation factor must be between 0.1 and 10.0.");

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
            {
                salinityField[i] = 35.0 * (1.0 - Math.Exp(-positions[i] / (estuaryLength / 2.0))); // Adjusted for broader ETM
                shearRate[i] = 0.0;
                flocSizes[i] = grainSizes[2];
            }
            settlingVelocities = CalculateSettlingVelocities();
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
            double tidalVelocity = 1.5 * tidalAmplitude * Math.Sin(2 * Math.PI * time / tidalPeriod); // Amplified tidal effect
            double riverVelocity = riverInflow / (estuaryDepth * (estuaryLength / dx)) * 0.5; // Reduced base flow
            for (int i = 0; i < gridPoints; i++)
            {
                double x = positions[i];
                double distanceFactor = Math.Exp(-x / estuaryLength);
                velocityField[i] = riverVelocity + tidalVelocity * distanceFactor;
            }
        }

        private void UpdateSimulation()
        {
            UpdateVelocityField();
            double dx = estuaryLength / (gridPoints - 1);
            double[] currentFluxes = new double[sedimentTypes.Length * gridPoints];
            double waterDensity = 1000.0; // kg/m³
            double waterViscosity = 1.0e-6; // m²/s
            double dragCoefficient = 0.0025;

            // Update shear rate and flocculation
            double[] dynamicSettlingVelocities = (double[])settlingVelocities.Clone();
            for (int i = 0; i < gridPoints; i++)
            {
                double shearStress = waterDensity * dragCoefficient * velocityField[i] * velocityField[i];
                shearRate[i] = Math.Sqrt(shearStress / (waterDensity * waterViscosity));
                double salinityFactor = Math.Exp(-Math.Pow((salinityField[i] - optimalSalinity) / salinityScale, 2));
                double shearFactor = optimalShearRate / (shearRate[i] + optimalShearRate) * Math.Exp(-shearRate[i] / breakupShearRate);
                double concentrationFactor = Math.Pow(Math.Max(sedimentConcentrations[2][i], 1e-6) / referenceConcentration, concentrationExponent);
                double flocFactor = Math.Max(flocculationCoefficient * concentrationFactor * salinityFactor * shearFactor, minFlocFactor);
                flocFactor = Math.Min(flocFactor, 100.0);
                dynamicSettlingVelocities[2] = settlingVelocities[2] * flocFactor * 2.0; // Amplify clay settling
                dynamicSettlingVelocities[2] = Math.Min(dynamicSettlingVelocities[2], 0.01);
                flocSizes[i] = flocSizeBase * Math.Pow(flocFactor, 1.0 / 3.0);
                flocSizes[i] = Math.Min(flocSizes[i], 0.0005);
            }

            // Update sediment concentrations
            for (int f = 0; f < sedimentTypes.Length; f++)
            {
                double[] newConcentrations = new double[gridPoints];
                double settlingVelocity = f == 2 ? dynamicSettlingVelocities[2] : settlingVelocities[f];
                for (int i = 1; i < gridPoints - 1; i++)
                {
                    double advection = velocityField[i] >= 0
                        ? -velocityField[i] * (sedimentConcentrations[f][i] - sedimentConcentrations[f][i - 1]) / dx
                        : -velocityField[i] * (sedimentConcentrations[f][i + 1] - sedimentConcentrations[f][i]) / dx;
                    double diffusion = diffusionCoefficients[f] * (sedimentConcentrations[f][i + 1] - 2 * sedimentConcentrations[f][i] + sedimentConcentrations[f][i - 1]) / (dx * dx);
                    double settling = -settlingVelocity * sedimentConcentrations[f][i];
                    newConcentrations[i] = sedimentConcentrations[f][i] + timeStep * (advection + diffusion + settling);
                    if (newConcentrations[i] < 0) newConcentrations[i] = 0;

                    double depositionRate = settlingVelocity * sedimentConcentrations[f][i] * (f == 2 && salinityField[i] >= 5.0 && salinityField[i] <= 15.0 ? 3.0 : 1.0); // Stronger ETM deposition
                    double shearStress = waterDensity * dragCoefficient * velocityField[i] * velocityField[i];
                    double erosionFactor = salinityField[i] >= 5.0 && salinityField[i] <= 15.0 && Math.Abs(velocityField[i]) < 0.5 ? 0.1 : 0.5; // Lower erosion during slack tides
                    double erosionRate = shearStress > criticalShearStresses[f] ? erosionFactor * 0.003 * (shearStress - criticalShearStresses[f]) : 0.0;
                    currentFluxes[f * gridPoints + i] = depositionRate - erosionRate;
                }
                newConcentrations[0] = sedimentConcentrations[f][0];
                newConcentrations[gridPoints - 1] = 0.0;
                sedimentConcentrations[f] = newConcentrations;
            }

            // Update salinity
            double[] newSalinity = new double[gridPoints];
            double tidalPhase = 2 * Math.PI * time / tidalPeriod;
            double intrusionFactor = 0.5 * (1 + Math.Cos(tidalPhase));
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
            newSalinity[0] = 0.0;
            newSalinity[gridPoints - 1] = 35.0 * (0.8 + 0.2 * intrusionFactor);
            salinityField = newSalinity;

            sedimentFluxes.Add(currentFluxes);
            if (sedimentFluxes.Count > maxTimeSteps)
                sedimentFluxes.RemoveAt(0);

            time += timeStep;
            visualizationPanel.Invalidate();
            if (showFluxCheckBox.Checked)
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
            double minVelocity = velocityField[0];
            double maxSalinity = salinityField[0];
            double minSalinity = salinityField[0];
            double maxETMConc = 0.0;
            double maxFlocSize = flocSizes[0];
            double maxFlocSettling = settlingVelocities[2];
            double avgShearRateETM = 0.0;
            int etmCount = 0;
            int etmIndex = -1;
            double[] maxFluxes = new double[sedimentTypes.Length];
            if (sedimentFluxes.Count > 0)
            {
                var latestFluxes = sedimentFluxes[sedimentFluxes.Count - 1];
                for (int f = 0; f < sedimentTypes.Length; f++)
                {
                    for (int i = 0; i < gridPoints; i++)
                    {
                        maxFluxes[f] = Math.Max(maxFluxes[f], Math.Abs(latestFluxes[f * gridPoints + i]));
                    }
                }
            }
            for (int i = 0; i < gridPoints; i++)
            {
                maxVelocity = Math.Max(maxVelocity, velocityField[i]);
                minVelocity = Math.Min(minVelocity, velocityField[i]);
                maxSalinity = Math.Max(maxSalinity, salinityField[i]);
                minSalinity = Math.Min(minSalinity, salinityField[i]);
                maxFlocSize = Math.Max(maxFlocSize, flocSizes[i]);
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
                    double salinityFactor = Math.Exp(-Math.Pow((salinityField[i] - optimalSalinity) / salinityScale, 2));
                    double shearFactor = optimalShearRate / (shearRate[i] + optimalShearRate) * Math.Exp(-shearRate[i] / breakupShearRate);
                    double concentrationFactor = Math.Pow(Math.Max(sedimentConcentrations[2][i], 1e-6) / referenceConcentration, concentrationExponent);
                    double flocFactor = Math.Max(flocculationCoefficient * concentrationFactor * salinityFactor * shearFactor, minFlocFactor);
                    flocFactor = Math.Min(flocFactor, 100.0);
                    double flocSettling = settlingVelocities[2] * flocFactor * 2.0;
                    maxFlocSettling = Math.Max(maxFlocSettling, flocSettling);
                    avgShearRateETM += shearRate[i];
                    etmCount++;
                }
            }
            avgShearRateETM = etmCount > 0 ? avgShearRateETM / etmCount : 0.0;
            maxFlocSettling = Math.Min(maxFlocSettling, 0.01);
            string output = $"Time: {time:F2}s | Tidal Phase: {phaseStr} | Vel (range): [{minVelocity:F6}, {maxVelocity:F6}] m/s | Salinity (range): [{minSalinity:F2}, {maxSalinity:F2}] ppt | ETM Conc: {maxETMConc:F6} kg/m³ @ x={(etmIndex >= 0 ? positions[etmIndex] : 0):F0}m | Max Floc Size: {(maxFlocSize * 1e6):F2} µm | Max Floc Settling: {maxFlocSettling:F6} m/s | Avg ETM Shear: {avgShearRateETM:F2} s^-1 | Fluxes: [Sand: {maxFluxes[0]:F6}, Silt: {maxFluxes[1]:F6}, Clay: {maxFluxes[2]:F6}] kg/m²s\r\n";
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
                double maxVelocity = 2.0; // Increased to accommodate tidal variations
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
            if (salinityRange < 1e-6) salinityRange = 1.0;
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

            if (sedimentFluxes.Count == 0)
            {
                g.DrawString("No flux data available", new Font("Verdana", 11), Brushes.Black, 50, 50);
                return;
            }

            int panelWidth = hovmoellerPanel.Width;
            int panelHeight = hovmoellerPanel.Height;
            int selectedSediment = sedimentSelector.SelectedIndex == 3 ? 2 : sedimentSelector.SelectedIndex; // Default to clay for "All"
            var latestFluxes = sedimentFluxes[sedimentFluxes.Count - 1];

            // Find flux range
            double maxFlux = 1e-6;
            double minFlux = -1e-6;
            for (int i = 0; i < gridPoints; i++)
            {
                double flux = latestFluxes[selectedSediment * gridPoints + i];
                maxFlux = Math.Max(maxFlux, flux);
                minFlux = Math.Min(minFlux, flux);
            }
            double fluxRange = Math.Max(maxFlux, -minFlux);
            if (fluxRange < 1e-10) fluxRange = 1e-6;

            // Define plot area with margins
            float plotX = 50;
            float plotY = 50;
            float plotWidth = panelWidth - 100;
            float plotHeight = panelHeight - 100;
            g.DrawRectangle(Pens.Black, plotX, plotY, plotWidth, plotHeight);

            // Highlight ETM zone (5–15 ppt)
            int etmStart = -1, etmEnd = -1;
            for (int i = 0; i < gridPoints; i++)
            {
                if (salinityField[i] >= 5.0 && salinityField[i] <= 15.0)
                {
                    if (etmStart == -1) etmStart = i;
                    etmEnd = i;
                }
            }
            if (etmStart != -1)
            {
                float xStart = plotX + etmStart * plotWidth / (float)(gridPoints - 1);
                float xEnd = plotX + etmEnd * plotWidth / (float)(gridPoints - 1);
                g.FillRectangle(new SolidBrush(Color.FromArgb(50, 0, 255, 0)), xStart, plotY, xEnd - xStart, plotHeight);
            }

            // Draw zero line
            float zeroY = plotY + plotHeight / 2;
            g.DrawLine(Pens.Black, plotX, zeroY, plotX + plotWidth, zeroY);

            // Plot flux curve
            using (Pen pen = new Pen(selectedSediment == 0 ? Color.Red : selectedSediment == 1 ? Color.Blue : Color.Green, 2))
            {
                for (int i = 0; i < gridPoints - 1; i++)
                {
                    float x1 = plotX + i * plotWidth / (float)(gridPoints - 1);
                    float x2 = plotX + (i + 1) * plotWidth / (float)(gridPoints - 1);
                    double flux1 = latestFluxes[selectedSediment * gridPoints + i];
                    double flux2 = latestFluxes[selectedSediment * gridPoints + i + 1];
                    float y1 = zeroY - (float)(flux1 / fluxRange * (plotHeight / 2));
                    float y2 = zeroY - (float)(flux2 / fluxRange * (plotHeight / 2));
                    g.DrawLine(pen, x1, y1, x2, y2);
                }
            }

            // Draw axes and labels
            g.DrawString("Distance (m)", new Font("Verdana", 11), Brushes.Black, plotX + plotWidth / 2 - 50, plotY + plotHeight + 10);
            g.DrawString($"{sedimentTypes[selectedSediment]} Flux (kg/m²s)", new Font("Verdana", 11), Brushes.Black, plotX - 50, plotY + plotHeight / 2 - 20, new StringFormat { Alignment = StringAlignment.Far, FormatFlags = StringFormatFlags.DirectionVertical });
            g.DrawString($"Time: {time:F2} s | Phase: {(2 * Math.PI * time / tidalPeriod % (2 * Math.PI) < Math.PI ? "Flood" : "Ebb")}", new Font("Verdana", 11), Brushes.Black, plotX, plotY - 20);
            g.DrawString($"Flux Range: [{minFlux:F6}, {maxFlux:F6}] kg/m²s", new Font("Verdana", 11), Brushes.Black, plotX + plotWidth - 200, plotY - 20);
            if (etmStart != -1)
            {
                float xStart = plotX + etmStart * plotWidth / (float)(gridPoints - 1);
                float xEnd = plotX + etmEnd * plotWidth / (float)(gridPoints - 1);
                g.DrawString("ETM Zone", new Font("Verdana", 10), Brushes.Black, xStart + (xEnd - xStart) / 2 - 30, plotY + 10);
            }

            // Legend
            g.FillRectangle(new SolidBrush(Color.FromArgb(50, 0, 255, 0)), plotX + plotWidth - 70, plotY + 20, 15, 15);
            g.DrawString("ETM Zone (5-15 ppt)", new Font("Verdana", 11), Brushes.Black, plotX + plotWidth - 50, plotY + 20);
            using (Pen pen = new Pen(selectedSediment == 0 ? Color.Red : selectedSediment == 1 ? Color.Blue : Color.Green, 2))
            {
                g.DrawLine(pen, plotX + plotWidth - 70, plotY + 40, plotX + plotWidth - 50, plotY + 40);
            }
            g.DrawString(sedimentTypes[selectedSediment], new Font("Verdana", 11), Brushes.Black, plotX + plotWidth - 50, plotY + 35);
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
            double maxFlocSize = flocSizes[0];
            for (int i = 0; i < gridPoints; i++)
            {
                minSalinity = Math.Min(minSalinity, salinityField[i]);
                maxSalinity = Math.Max(maxSalinity, salinityField[i]);
                maxFlocSize = Math.Max(maxFlocSize, flocSizes[i]);
            }
            double salinityRange = maxSalinity - minSalinity;
            if (salinityRange < 1e-6) salinityRange = 35.0;
            if (maxFlocSize < 1e-6) maxFlocSize = 0.0005;

            if (showAllSediments)
            {
                // Correlation matrix for all sediments
                int gridSize = panelWidth / 3;
                double maxETMConc = 0.0;
                double maxETMFloc = 0.0;
                for (int i = 0; i < gridPoints; i++)
                {
                    if (salinityField[i] >= 5.0 && salinityField[i] <= 15.0)
                    {
                        for (int f = 0; f < sedimentTypes.Length; f++)
                            maxETMConc = Math.Max(maxETMConc, sedimentConcentrations[f][i]);
                        maxETMFloc = Math.Max(maxETMFloc, flocSizes[i]);
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
                        int size = f == 2 ? (int)(5 + 15 * flocSizes[i] / maxFlocSize) : (int)(5 + 10 * sedimentConcentrations[f][i] / maxConcentration);
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

                // ETM concentration and floc size lines
                float etmY = panelHeight - (float)(maxETMConc / maxConcentration * (panelHeight - 100)) + 50;
                g.DrawLine(new Pen(Color.Black, 1) { DashStyle = System.Drawing.Drawing2D.DashStyle.Dash }, 50, etmY, panelWidth - 50, etmY);
                g.DrawString($"Max ETM Conc: {maxETMConc:F6} kg/m³", new Font("Verdana", 10), Brushes.Black, panelWidth - 200, etmY - 20);
                if (maxETMFloc > 0)
                {
                    float flocY = panelHeight - (float)(maxETMFloc / maxFlocSize * (panelHeight - 100)) + 50;
                    g.DrawLine(new Pen(Color.Gray, 1) { DashStyle = System.Drawing.Drawing2D.DashStyle.Dot }, 50, flocY, panelWidth - 50, flocY);
                    g.DrawString($"Max ETM Floc: {(maxETMFloc * 1e6):F2} µm", new Font("Verdana", 10), Brushes.Black, panelWidth - 200, flocY - 40);
                }

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
                g.FillEllipse(new SolidBrush(Color.Green), panelWidth - 70, 100, 15, 15);
                g.DrawString("Clay Floc Size", new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 100);
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
                double maxETMFloc = 0.0;
                for (int i = 0; i < gridPoints; i++)
                {
                    if (salinityField[i] >= 5.0 && salinityField[i] <= 15.0)
                    {
                        maxETMConc = Math.Max(maxETMConc, sedimentConcentrations[selectedSediment][i]);
                        if (selectedSediment == 2)
                            maxETMFloc = Math.Max(maxETMFloc, flocSizes[i]);
                    }
                    float x = plotX + (float)((salinityField[i] - minSalinity) / salinityRange * plotWidth);
                    float y = plotY + plotHeight - (float)(sedimentConcentrations[selectedSediment][i] / maxConcentration * plotHeight);
                    int size = selectedSediment == 2 ? (int)(5 + 15 * flocSizes[i] / maxFlocSize) : (int)(5 + 10 * sedimentConcentrations[selectedSediment][i] / maxConcentration);
                    int alpha = salinityField[i] >= 5.0 && salinityField[i] <= 15.0 ? 255 : 128;
                    using (Brush brush = new SolidBrush(Color.FromArgb(alpha, color)))
                    {
                        g.FillEllipse(brush, x - size / 2, y - size / 2, size, size);
                    }
                }

                // ETM concentration and floc size lines
                float etmY = plotY + plotHeight - (float)(maxETMConc / maxConcentration * plotHeight);
                g.DrawLine(new Pen(Color.Black, 1) { DashStyle = System.Drawing.Drawing2D.DashStyle.Dash }, plotX, etmY, plotX + plotWidth, etmY);
                g.DrawString($"Max ETM Conc: {maxETMConc:F6} kg/m³", new Font("Verdana", 10), Brushes.Black, plotX + plotWidth - 200, etmY - 20);
                if (selectedSediment == 2 && maxETMFloc > 0)
                {
                    float flocY = plotY + plotHeight - (float)(maxETMFloc / maxFlocSize * plotHeight);
                    g.DrawLine(new Pen(Color.Gray, 1) { DashStyle = System.Drawing.Drawing2D.DashStyle.Dot }, plotX, flocY, plotX + plotWidth, flocY);
                    g.DrawString($"Max ETM Floc: {(maxETMFloc * 1e6):F2} µm", new Font("Verdana", 10), Brushes.Black, plotX + plotWidth - 200, flocY - 40);
                }

                g.DrawString($"{sedimentTypes[selectedSediment]} vs. Salinity", new Font("Verdana", 11), Brushes.Black, plotX + 5, plotY - 20);
                g.DrawString($"Salinity (ppt) [{minSalinity:F2}, {maxSalinity:F2}]", new Font("Verdana", 11), Brushes.Black, plotX, plotY + plotHeight + 5);
                g.DrawString("Concentration (kg/m³)", new Font("Verdana", 11), Brushes.Black, plotX - 50, plotY + plotHeight / 2);

                // Legend
                g.FillRectangle(new SolidBrush(Color.FromArgb(50, 0, 255, 0)), panelWidth - 70, 20, 15, 15);
                g.DrawString("ETM Zone (5-15 ppt)", new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 20);
                g.FillEllipse(new SolidBrush(color), panelWidth - 70, 40, 10, 10);
                g.DrawString(sedimentTypes[selectedSediment], new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 40);
                if (selectedSediment == 2)
                {
                    g.FillEllipse(new SolidBrush(Color.Green), panelWidth - 70, 60, 15, 15);
                    g.DrawString("Clay Floc Size", new Font("Verdana", 11), Brushes.Black, panelWidth - 50, 60);
                }
            }
        }

        private void UpdateButtonStates()
        {
            startSimButton.Enabled = !isRunning;
            pauseSimButton.Enabled = isRunning;
            resetSimButton.Enabled = true;
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
            shearRate = new double[gridPoints];
            flocSizes = new double[gridPoints];
            sedimentFluxes = new List<double[]>();
            for (int i = 0; i < gridPoints; i++)
            {
                positions[i] = i * estuaryLength / (gridPoints - 1);
                velocityField[i] = riverInflow / (estuaryDepth * (estuaryLength / (estuaryLength / (gridPoints - 1))));
                salinityField[i] = 35.0 * (1.0 - Math.Exp(-positions[i] / (estuaryLength / 2.0))); // Adjusted for broader ETM
                shearRate[i] = 0.0;
                flocSizes[i] = grainSizes[2]; // Initial clay floc size = 2 µm
            }

            simulationTimer = new Timer
            {
                Interval = 100
            };
            simulationTimer.Tick += (s, e) => UpdateSimulation();
            UpdateButtonStates();
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