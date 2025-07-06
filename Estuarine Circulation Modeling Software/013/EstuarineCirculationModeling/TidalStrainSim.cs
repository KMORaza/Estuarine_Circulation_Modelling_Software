using System;
using System.Drawing;
using System.Windows.Forms;

namespace EstuarineCirculationModeling
{
    public class TidalStrainSim
    {
        private Form tidalStrainForm;
        private PictureBox visualizationBox;
        private PictureBox timeSeriesBox;
        private Label timeLabel;
        private TextBox timeTextBox;
        private Label strainRateLabel;
        private TextBox strainRateTextBox;
        private Label richardsonLabel;
        private TextBox richardsonTextBox;
        private Label viscosityLabel;
        private TextBox viscosityTextBox;
        private Label thermalDiffusivityLabel;
        private TextBox thermalDiffusivityTextBox;
        private Label sedimentDiffusivityLabel;
        private TextBox sedimentDiffusivityTextBox;
        private Label settlingVelocityLabel;
        private TextBox settlingVelocityTextBox;
        private Label coriolisLabel;
        private TextBox coriolisTextBox;
        private Label turbulenceLabel;
        private TextBox turbulenceTextBox;
        private Label tidalStrainLabel;
        private TextBox tidalStrainTextBox;
        private Label densityRatioLabel;
        private TextBox densityRatioTextBox;
        private ComboBox visualizationModeComboBox;
        private Button startButton;
        private Button pauseButton;
        private Button resetButton;
        private Timer simulationTimer;
        private double currentTime;
        private double tidalPeriod;
        private double tidalAmplitude;
        private double estuaryLength;
        private double estuaryDepth;
        private double[,] salinityProfile;
        private double[,] temperatureProfile;
        private double[,] sedimentProfile; // Suspended sediment concentration (g/m³)
        private double[] bedloadProfile; // Bedload mass (kg/m²)
        private double[,] uVelocityProfile; // Horizontal velocity
        private double[,] wVelocityProfile; // Vertical velocity
        private double kinematicViscosity;
        private double thermalDiffusivity;
        private double sedimentDiffusivity;
        private double settlingVelocity;
        private double coriolisParameter;
        private double turbulenceIntensity;
        private bool isSimulationRunning;
        private const int gridX = 400; // Finer grid for DNS
        private const int gridZ = 100; // Vertical grid for depth
        private readonly Random rand = new Random();
        private float zoomLevel = 1.0f;
        private readonly double[] richardsonHistory = new double[1000];
        private readonly double[] densityRatioHistory = new double[1000];
        private readonly double[] sedimentHistory = new double[1000];
        private int historyIndex = 0;
        private const double alpha = 2e-4; // Thermal expansion coefficient (1/°C)
        private const double beta = 7.6e-4; // Saline contraction coefficient (1/PSU)
        private const double rho0 = 1000.0; // Reference density (kg/m³)
        private const double rhoS = 2650.0; // Sediment density (kg/m³)
        private const double gamma = (rhoS - rho0) / (rho0 * rhoS); // Sediment density coefficient
        private const double T0 = 15.0; // Reference temperature (°C)
        private const double S0 = 35.0; // Reference salinity (PSU)
        private const double C0 = 100.0; // Reference sediment concentration (g/m³)
        private const double g = 9.81; // Gravity (m/s²)
        private const double d = 0.001; // Grain size (m)
        private const double thetaC = 0.047; // Critical Shields parameter
        private const double alphaE = 0.0001; // Erosion coefficient

        public TidalStrainSim(double estuaryLength, double estuaryDepth, double tidalPeriod, double tidalAmplitude)
        {
            this.estuaryLength = estuaryLength;
            this.estuaryDepth = estuaryDepth;
            this.tidalPeriod = tidalPeriod;
            this.tidalAmplitude = tidalAmplitude;
            this.currentTime = 0.0;
            this.kinematicViscosity = 1e-6; // m²/s (water viscosity)
            this.thermalDiffusivity = 1.4e-7; // m²/s (thermal diffusivity)
            this.sedimentDiffusivity = 1e-6; // m²/s (sediment diffusivity)
            this.settlingVelocity = -0.001; // m/s (negative for sinking)
            this.coriolisParameter = 1e-4; // s^-1 (mid-latitude)
            this.turbulenceIntensity = 0.01; // Turbulence scale
            this.isSimulationRunning = false;
            this.salinityProfile = new double[gridX, gridZ];
            this.temperatureProfile = new double[gridX, gridZ];
            this.sedimentProfile = new double[gridX, gridZ];
            this.bedloadProfile = new double[gridX];
            this.uVelocityProfile = new double[gridX, gridZ];
            this.wVelocityProfile = new double[gridX, gridZ];
            InitializeProfiles();
            InitializeForm();
            InitializeSimulationTimer();
        }

        private void InitializeProfiles()
        {
            // Initialize salinity, temperature, sediment, and bedload
            for (int i = 0; i < gridX; i++)
            {
                for (int j = 0; j < gridZ; j++)
                {
                    double xFraction = (double)i / (gridX - 1);
                    double zFraction = (double)j / (gridZ - 1);
                    salinityProfile[i, j] = 35.0 * xFraction * (1.0 - 0.5 * zFraction); // Stratified salinity
                    temperatureProfile[i, j] = 15.0 - 5.0 * zFraction + 2.0 * xFraction; // 15°C at surface, 10°C at bottom
                    sedimentProfile[i, j] = 100.0 * (1.0 - zFraction); // Higher concentration near bottom
                }
                bedloadProfile[i] = 0.1; // Initial uniform bedload (kg/m²)
            }
        }

        private void InitializeForm()
        {
            tidalStrainForm = new Form
            {
                Text = "Tidal Straining Simulation (DNS)",
                Size = new Size(1100, 800),
                FormBorderStyle = FormBorderStyle.FixedDialog,
                MaximizeBox = false,
                Font = new Font("Consolas", 9F)
            };

            // Control Panel
            Panel controlPanel = new Panel
            {
                Location = new Point(10, 10),
                Size = new Size(200, 740),
                BorderStyle = BorderStyle.FixedSingle,
                AutoScroll = true
            };

            // Time Label and TextBox
            timeLabel = new Label
            {
                Location = new Point(10, 20),
                Size = new Size(80, 15),
                Text = "Time (s):",
                AutoSize = true
            };
            timeTextBox = new TextBox
            {
                Location = new Point(10, 40),
                Size = new Size(150, 22),
                ReadOnly = true,
                Text = "0.0"
            };

            // Strain Rate Label and TextBox
            strainRateLabel = new Label
            {
                Location = new Point(10, 70),
                Size = new Size(80, 15),
                Text = "Strain Rate (1/s):",
                AutoSize = true
            };
            strainRateTextBox = new TextBox
            {
                Location = new Point(10, 90),
                Size = new Size(150, 22),
                ReadOnly = true,
                Text = "0.0"
            };

            // Richardson Number Label and TextBox
            richardsonLabel = new Label
            {
                Location = new Point(10, 120),
                Size = new Size(100, 15),
                Text = "Avg Richardson #:",
                AutoSize = true
            };
            richardsonTextBox = new TextBox
            {
                Location = new Point(10, 140),
                Size = new Size(150, 22),
                ReadOnly = true,
                Text = "0.0"
            };

            // Viscosity Label and TextBox
            viscosityLabel = new Label
            {
                Location = new Point(10, 170),
                Size = new Size(100, 15),
                Text = "Viscosity (m²/s):",
                AutoSize = true
            };
            viscosityTextBox = new TextBox
            {
                Location = new Point(10, 190),
                Size = new Size(150, 22),
                Text = kinematicViscosity.ToString("E2")
            };
            viscosityTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(viscosityTextBox.Text, out double newViscosity))
                {
                    kinematicViscosity = Math.Max(1e-7, Math.Min(1e-5, newViscosity));
                    sedimentDiffusivity = kinematicViscosity; // Sync with viscosity
                    sedimentDiffusivityTextBox.Text = sedimentDiffusivity.ToString("E2");
                }
            };

            // Thermal Diffusivity Label and TextBox
            thermalDiffusivityLabel = new Label
            {
                Location = new Point(10, 220),
                Size = new Size(100, 15),
                Text = "Therm. Diff. (m²/s):",
                AutoSize = true
            };
            thermalDiffusivityTextBox = new TextBox
            {
                Location = new Point(10, 240),
                Size = new Size(150, 22),
                Text = thermalDiffusivity.ToString("E2")
            };
            thermalDiffusivityTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(thermalDiffusivityTextBox.Text, out double newDiffusivity))
                {
                    thermalDiffusivity = Math.Max(1e-8, Math.Min(1e-6, newDiffusivity));
                }
            };

            // Sediment Diffusivity Label and TextBox
            sedimentDiffusivityLabel = new Label
            {
                Location = new Point(10, 270),
                Size = new Size(100, 15),
                Text = "Sed. Diff. (m²/s):",
                AutoSize = true
            };
            sedimentDiffusivityTextBox = new TextBox
            {
                Location = new Point(10, 290),
                Size = new Size(150, 22),
                Text = sedimentDiffusivity.ToString("E2")
            };
            sedimentDiffusivityTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(sedimentDiffusivityTextBox.Text, out double newDiffusivity))
                {
                    sedimentDiffusivity = Math.Max(1e-7, Math.Min(1e-5, newDiffusivity));
                }
            };

            // Settling Velocity Label and TextBox
            settlingVelocityLabel = new Label
            {
                Location = new Point(10, 320),
                Size = new Size(100, 15),
                Text = "Settle Vel. (m/s):",
                AutoSize = true
            };
            settlingVelocityTextBox = new TextBox
            {
                Location = new Point(10, 340),
                Size = new Size(150, 22),
                Text = settlingVelocity.ToString("E2")
            };
            settlingVelocityTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(settlingVelocityTextBox.Text, out double newVelocity))
                {
                    settlingVelocity = Math.Min(0.0, Math.Max(-0.01, newVelocity)); // Negative, capped
                }
            };

            // Coriolis Label and TextBox
            coriolisLabel = new Label
            {
                Location = new Point(10, 370),
                Size = new Size(100, 15),
                Text = "Coriolis (1/s):",
                AutoSize = true
            };
            coriolisTextBox = new TextBox
            {
                Location = new Point(10, 390),
                Size = new Size(150, 22),
                Text = coriolisParameter.ToString("E2")
            };
            coriolisTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(coriolisTextBox.Text, out double newCoriolis))
                {
                    coriolisParameter = Math.Max(0.0, Math.Min(1e-3, newCoriolis));
                }
            };

            // Turbulence Intensity Label and TextBox
            turbulenceLabel = new Label
            {
                Location = new Point(10, 420),
                Size = new Size(100, 15),
                Text = "Turbulence Intensity:",
                AutoSize = true
            };
            turbulenceTextBox = new TextBox
            {
                Location = new Point(10, 440),
                Size = new Size(150, 22),
                Text = turbulenceIntensity.ToString("F3")
            };
            turbulenceTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(turbulenceTextBox.Text, out double newTurbulence))
                {
                    turbulenceIntensity = Math.Max(0.0, Math.Min(0.1, newTurbulence));
                }
            };

            // Tidal Strain Term Label and TextBox
            tidalStrainLabel = new Label
            {
                Location = new Point(10, 470),
                Size = new Size(100, 15),
                Text = "Tidal Strain (PSU/s):",
                AutoSize = true
            };
            tidalStrainTextBox = new TextBox
            {
                Location = new Point(10, 490),
                Size = new Size(150, 22),
                ReadOnly = true,
                Text = "0.0"
            };

            // Density Ratio Label and TextBox
            densityRatioLabel = new Label
            {
                Location = new Point(10, 520),
                Size = new Size(100, 15),
                Text = "Density Ratio:",
                AutoSize = true
            };
            densityRatioTextBox = new TextBox
            {
                Location = new Point(10, 540),
                Size = new Size(150, 22),
                ReadOnly = true,
                Text = "0.0"
            };

            // Visualization Mode ComboBox
            visualizationModeComboBox = new ComboBox
            {
                Location = new Point(10, 570),
                Size = new Size(150, 22),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            visualizationModeComboBox.Items.AddRange(new[] { "Salinity", "Vorticity", "Density", "Turbidity", "Bedload" });
            visualizationModeComboBox.SelectedIndex = 0;
            visualizationModeComboBox.SelectedIndexChanged += (s, e) => visualizationBox.Invalidate();

            // Buttons
            startButton = new Button
            {
                Location = new Point(10, 600),
                Size = new Size(150, 25),
                Text = "Start",
                FlatStyle = FlatStyle.Flat
            };
            startButton.Click += StartButton_Click;

            pauseButton = new Button
            {
                Location = new Point(10, 635),
                Size = new Size(150, 25),
                Text = "Pause",
                FlatStyle = FlatStyle.Flat,
                Enabled = false
            };
            pauseButton.Click += PauseButton_Click;

            resetButton = new Button
            {
                Location = new Point(10, 670),
                Size = new Size(150, 25),
                Text = "Reset",
                FlatStyle = FlatStyle.Flat
            };
            resetButton.Click += ResetButton_Click;

            // Visualization Area
            visualizationBox = new PictureBox
            {
                Location = new Point(220, 10),
                Size = new Size(650, 450),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            visualizationBox.Paint += VisualizationBox_Paint;
            visualizationBox.MouseWheel += (s, e) =>
            {
                zoomLevel = Math.Max(0.5f, Math.Min(2.0f, zoomLevel + e.Delta * 0.001f));
                visualizationBox.Invalidate();
            };

            // Time Series Plot
            timeSeriesBox = new PictureBox
            {
                Location = new Point(220, 470),
                Size = new Size(650, 300),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            timeSeriesBox.Paint += TimeSeriesBox_Paint;

            // Add controls to panel
            controlPanel.Controls.Add(timeLabel);
            controlPanel.Controls.Add(timeTextBox);
            controlPanel.Controls.Add(strainRateLabel);
            controlPanel.Controls.Add(strainRateTextBox);
            controlPanel.Controls.Add(richardsonLabel);
            controlPanel.Controls.Add(richardsonTextBox);
            controlPanel.Controls.Add(viscosityLabel);
            controlPanel.Controls.Add(viscosityTextBox);
            controlPanel.Controls.Add(thermalDiffusivityLabel);
            controlPanel.Controls.Add(thermalDiffusivityTextBox);
            controlPanel.Controls.Add(sedimentDiffusivityLabel);
            controlPanel.Controls.Add(sedimentDiffusivityTextBox);
            controlPanel.Controls.Add(settlingVelocityLabel);
            controlPanel.Controls.Add(settlingVelocityTextBox);
            controlPanel.Controls.Add(coriolisLabel);
            controlPanel.Controls.Add(coriolisTextBox);
            controlPanel.Controls.Add(turbulenceLabel);
            controlPanel.Controls.Add(turbulenceTextBox);
            controlPanel.Controls.Add(tidalStrainLabel);
            controlPanel.Controls.Add(tidalStrainTextBox);
            controlPanel.Controls.Add(densityRatioLabel);
            controlPanel.Controls.Add(densityRatioTextBox);
            controlPanel.Controls.Add(visualizationModeComboBox);
            controlPanel.Controls.Add(startButton);
            controlPanel.Controls.Add(pauseButton);
            controlPanel.Controls.Add(resetButton);

            // Add controls to form
            tidalStrainForm.Controls.Add(controlPanel);
            tidalStrainForm.Controls.Add(visualizationBox);
            tidalStrainForm.Controls.Add(timeSeriesBox);

            // Handle form closing to stop the timer
            tidalStrainForm.FormClosing += (s, e) =>
            {
                simulationTimer.Stop();
                isSimulationRunning = false;
                UpdateButtonStates();
            };
        }

        private void InitializeSimulationTimer()
        {
            simulationTimer = new Timer();
            simulationTimer.Interval = 10; // 10ms for finer time steps
            simulationTimer.Tick += (s, e) =>
            {
                UpdateSimulation();
                visualizationBox.Invalidate();
                timeSeriesBox.Invalidate();
            };
        }

        private double ComputeDensity(int i, int j)
        {
            return rho0 * (1.0 - alpha * (temperatureProfile[i, j] - T0) + beta * (salinityProfile[i, j] - S0) + gamma * sedimentProfile[i, j]);
        }

        private void UpdateSimulation()
        {
            currentTime += simulationTimer.Interval / 1000.0; // Convert ms to s
            double tidalPhase = 2 * Math.PI * currentTime / tidalPeriod;
            double tidalVelocity = tidalAmplitude * Math.Cos(tidalPhase); // Base tidal velocity (m/s)

            double dx = estuaryLength / (gridX - 1);
            double dz = estuaryDepth / (gridZ - 1);
            double dt = simulationTimer.Interval / 1000.0; // 0.01s

            // Update velocity profiles (Navier-Stokes: advection + viscosity + Coriolis)
            double[,] newUVelocityProfile = new double[gridX, gridZ];
            double[,] newWVelocityProfile = new double[gridX, gridZ];
            for (int i = 1; i < gridX - 1; i++)
            {
                for (int j = 1; j < gridZ - 1; j++)
                {
                    // Upwind advection for u: u * du/dx + w * du/dz
                    double dudx = uVelocityProfile[i, j] > 0 ?
                        (uVelocityProfile[i, j] - uVelocityProfile[i - 1, j]) / dx :
                        (uVelocityProfile[i + 1, j] - uVelocityProfile[i, j]) / dx;
                    double dudz = wVelocityProfile[i, j] > 0 ?
                        (uVelocityProfile[i, j] - uVelocityProfile[i, j - 1]) / dz :
                        (uVelocityProfile[i, j + 1] - uVelocityProfile[i, j]) / dz;
                    double uAdvection = uVelocityProfile[i, j] * dudx + wVelocityProfile[i, j] * dudz;

                    // Viscosity term: ν * (d²u/dx² + d²u/dz²)
                    double d2udx2 = (uVelocityProfile[i + 1, j] - 2 * uVelocityProfile[i, j] + uVelocityProfile[i - 1, j]) / (dx * dx);
                    double d2udz2 = (uVelocityProfile[i, j + 1] - 2 * uVelocityProfile[i, j] + uVelocityProfile[i, j - 1]) / (dz * dz);
                    double uViscosityTerm = kinematicViscosity * (d2udx2 + d2udz2);

                    // Coriolis term: -f * w
                    double uCoriolis = -coriolisParameter * wVelocityProfile[i, j];

                    // Tidal forcing
                    double uForcing = (tidalVelocity - uVelocityProfile[i, j]) / tidalPeriod;

                    newUVelocityProfile[i, j] = uVelocityProfile[i, j] + dt * (-uAdvection + uViscosityTerm + uCoriolis + uForcing);

                    // Upwind advection for w: u * dw/dx + w * dw/dz
                    double dwdx = uVelocityProfile[i, j] > 0 ?
                        (wVelocityProfile[i, j] - wVelocityProfile[i - 1, j]) / dx :
                        (wVelocityProfile[i + 1, j] - wVelocityProfile[i, j]) / dx;
                    double dwdz = wVelocityProfile[i, j] > 0 ?
                        (wVelocityProfile[i, j] - wVelocityProfile[i, j - 1]) / dz :
                        (wVelocityProfile[i, j + 1] - wVelocityProfile[i, j]) / dz;
                    double wAdvection = uVelocityProfile[i, j] * dwdx + wVelocityProfile[i, j] * dwdz;

                    // Viscosity term: ν * (d²w/dx² + d²w/dz²)
                    double d2wdx2 = (wVelocityProfile[i + 1, j] - 2 * wVelocityProfile[i, j] + wVelocityProfile[i - 1, j]) / (dx * dx);
                    double d2wdz2 = (wVelocityProfile[i, j + 1] - 2 * wVelocityProfile[i, j] + wVelocityProfile[i, j - 1]) / (dz * dz);
                    double wViscosityTerm = kinematicViscosity * (d2wdx2 + d2wdz2);

                    // Coriolis term: f * u
                    double wCoriolis = coriolisParameter * uVelocityProfile[i, j];

                    newWVelocityProfile[i, j] = wVelocityProfile[i, j] + dt * (-wAdvection + wViscosityTerm + wCoriolis);

                    // Add turbulence
                    newUVelocityProfile[i, j] += turbulenceIntensity * tidalAmplitude * (rand.NextDouble() - 0.5);
                    newWVelocityProfile[i, j] += turbulenceIntensity * tidalAmplitude * (rand.NextDouble() - 0.5);
                }
            }
            // Boundary conditions
            for (int i = 0; i < gridX; i++)
            {
                newUVelocityProfile[i, 0] = 0.0; // No-slip at bottom
                newUVelocityProfile[i, gridZ - 1] = tidalVelocity; // Ocean boundary
                newWVelocityProfile[i, 0] = 0.0; // No-slip at bottom
                newWVelocityProfile[i, gridZ - 1] = 0.0; // Free surface
            }
            for (int j = 0; j < gridZ; j++)
            {
                newUVelocityProfile[0, j] = 0.0; // River inflow
                newUVelocityProfile[gridX - 1, j] = tidalVelocity; // Ocean
                newWVelocityProfile[0, j] = 0.0; // River
                newWVelocityProfile[gridX - 1, j] = 0.0; // Ocean
            }
            uVelocityProfile = newUVelocityProfile;
            wVelocityProfile = newWVelocityProfile;

            // Update salinity profile (advection-diffusion)
            double[,] newSalinityProfile = new double[gridX, gridZ];
            double avgTidalStrain = 0.0;
            for (int i = 1; i < gridX - 1; i++)
            {
                for (int j = 1; j < gridZ - 1; j++)
                {
                    // Upwind advection: u * ds/dx + w * ds/dz
                    double dsdx = uVelocityProfile[i, j] > 0 ?
                        (salinityProfile[i, j] - salinityProfile[i - 1, j]) / dx :
                        (salinityProfile[i + 1, j] - salinityProfile[i, j]) / dx;
                    double dsdz = wVelocityProfile[i, j] > 0 ?
                        (salinityProfile[i, j] - salinityProfile[i, j - 1]) / dz :
                        (salinityProfile[i, j + 1] - salinityProfile[i, j]) / dz;
                    double advection = uVelocityProfile[i, j] * dsdx + wVelocityProfile[i, j] * dsdz;

                    // Diffusion term
                    double d2sdx2 = (salinityProfile[i + 1, j] - 2 * salinityProfile[i, j] + salinityProfile[i - 1, j]) / (dx * dx);
                    double d2sdz2 = (salinityProfile[i, j + 1] - 2 * salinityProfile[i, j] + salinityProfile[i, j - 1]) / (dz * dz);
                    double diffusion = kinematicViscosity * (d2sdx2 + d2sdz2);

                    newSalinityProfile[i, j] = salinityProfile[i, j] + dt * (-advection + diffusion);
                    avgTidalStrain += Math.Abs(uVelocityProfile[i, j] * dsdx);
                }
            }
            // Salinity boundary conditions
            for (int j = 0; j < gridZ; j++)
            {
                newSalinityProfile[0, j] = 0.0; // River
                newSalinityProfile[gridX - 1, j] = 35.0; // Ocean
            }
            for (int i = 0; i < gridX; i++)
            {
                newSalinityProfile[i, 0] = newSalinityProfile[i, 1]; // Bottom
                newSalinityProfile[i, gridZ - 1] = newSalinityProfile[i, gridZ - 2]; // Surface
            }
            salinityProfile = newSalinityProfile;

            // Update temperature profile (advection-diffusion)
            double[,] newTemperatureProfile = new double[gridX, gridZ];
            for (int i = 1; i < gridX - 1; i++)
            {
                for (int j = 1; j < gridZ - 1; j++)
                {
                    // Upwind advection: u * dT/dx + w * dT/dz
                    double dTdx = uVelocityProfile[i, j] > 0 ?
                        (temperatureProfile[i, j] - temperatureProfile[i - 1, j]) / dx :
                        (temperatureProfile[i + 1, j] - temperatureProfile[i, j]) / dx;
                    double dTdz = wVelocityProfile[i, j] > 0 ?
                        (temperatureProfile[i, j] - temperatureProfile[i, j - 1]) / dz :
                        (temperatureProfile[i, j + 1] - temperatureProfile[i, j]) / dz;
                    double advection = uVelocityProfile[i, j] * dTdx + wVelocityProfile[i, j] * dTdz;

                    // Diffusion term
                    double d2Tdx2 = (temperatureProfile[i + 1, j] - 2 * temperatureProfile[i, j] + temperatureProfile[i - 1, j]) / (dx * dx);
                    double d2Tdz2 = (temperatureProfile[i, j + 1] - 2 * temperatureProfile[i, j] + temperatureProfile[i, j - 1]) / (dz * dz);
                    double diffusion = thermalDiffusivity * (d2Tdx2 + d2Tdz2);

                    newTemperatureProfile[i, j] = temperatureProfile[i, j] + dt * (-advection + diffusion);
                }
            }
            // Temperature boundary conditions
            for (int j = 0; j < gridZ; j++)
            {
                newTemperatureProfile[0, j] = 12.0; // River (cooler inflow)
                newTemperatureProfile[gridX - 1, j] = 15.0; // Ocean (warmer)
            }
            for (int i = 0; i < gridX; i++)
            {
                newTemperatureProfile[i, 0] = newTemperatureProfile[i, 1]; // Bottom
                newTemperatureProfile[i, gridZ - 1] = 15.0; // Surface
            }
            temperatureProfile = newTemperatureProfile;

            // Update suspended sediment profile (advection-diffusion + settling)
            double[,] newSedimentProfile = new double[gridX, gridZ];
            double[] erosionRate = new double[gridX];
            double[] depositionRate = new double[gridX];
            for (int i = 1; i < gridX - 1; i++)
            {
                for (int j = 1; j < gridZ - 1; j++)
                {
                    // Upwind advection: u * dC/dx + (w + ws) * dC/dz
                    double dCdx = uVelocityProfile[i, j] > 0 ?
                        (sedimentProfile[i, j] - sedimentProfile[i - 1, j]) / dx :
                        (sedimentProfile[i + 1, j] - sedimentProfile[i, j]) / dx;
                    double dCdz = (wVelocityProfile[i, j] + settlingVelocity) > 0 ?
                        (sedimentProfile[i, j] - sedimentProfile[i, j - 1]) / dz :
                        (sedimentProfile[i, j + 1] - sedimentProfile[i, j]) / dz;
                    double advection = uVelocityProfile[i, j] * dCdx + (wVelocityProfile[i, j] + settlingVelocity) * dCdz;

                    // Diffusion term
                    double d2Cdx2 = (sedimentProfile[i + 1, j] - 2 * sedimentProfile[i, j] + sedimentProfile[i - 1, j]) / (dx * dx);
                    double d2Cdz2 = (sedimentProfile[i, j + 1] - 2 * sedimentProfile[i, j] + sedimentProfile[i, j - 1]) / (dz * dz);
                    double diffusion = sedimentDiffusivity * (d2Cdx2 + d2Cdz2);

                    newSedimentProfile[i, j] = Math.Max(0.0, sedimentProfile[i, j] + dt * (-advection + diffusion));
                }
            }
            // Compute bed shear stress, erosion, and deposition
            for (int i = 1; i < gridX - 1; i++)
            {
                double dudz = (uVelocityProfile[i, 1] - uVelocityProfile[i, 0]) / dz;
                double tauB = rho0 * kinematicViscosity * dudz;
                double uStar = Math.Sqrt(Math.Abs(tauB) / rho0);
                erosionRate[i] = alphaE * uStar * uStar;
                depositionRate[i] = -settlingVelocity * sedimentProfile[i, 1]; // Deposition at bottom
            }
            // Update bedload profile
            double[] newBedloadProfile = new double[gridX];
            for (int i = 1; i < gridX - 1; i++)
            {
                double dudz = (uVelocityProfile[i, 1] - uVelocityProfile[i, 0]) / dz;
                double tauB = rho0 * kinematicViscosity * dudz;
                double theta = tauB / ((rhoS - rho0) * g * d);
                double qb = theta > thetaC ? 8.0 * Math.Pow(theta - thetaC, 1.5) * Math.Sqrt((rhoS / rho0 - 1) * g * d * d * d) : 0.0;
                double dqb_dx = i < gridX - 2 ? (qb - (thetaC > thetaC ? 8.0 * Math.Pow(((rho0 * kinematicViscosity * (uVelocityProfile[i + 1, 1] - uVelocityProfile[i + 1, 0]) / dz) / ((rhoS - rho0) * g * d)) - thetaC, 1.5) * Math.Sqrt((rhoS / rho0 - 1) * g * d * d * d) : 0.0)) / dx : 0.0;
                newBedloadProfile[i] = Math.Max(0.0, bedloadProfile[i] + dt * (erosionRate[i] - depositionRate[i] - dqb_dx));
            }
            // Sediment boundary conditions
            for (int j = 0; j < gridZ; j++)
            {
                newSedimentProfile[0, j] = 0.0; // No sediment from river
                newSedimentProfile[gridX - 1, j] = 0.0; // No sediment from ocean
            }
            for (int i = 0; i < gridX; i++)
            {
                newSedimentProfile[i, 0] = newSedimentProfile[i, 1] + (erosionRate[i] - depositionRate[i]) * dt / dz; // Bottom exchange
                newSedimentProfile[i, gridZ - 1] = newSedimentProfile[i, gridZ - 2]; // Surface
                newBedloadProfile[i] = Math.Max(0.0, bedloadProfile[i]); // Ensure non-negative
            }
            sedimentProfile = newSedimentProfile;
            bedloadProfile = newBedloadProfile;

            avgTidalStrain /= (gridX - 2) * (gridZ - 2);
            tidalStrainTextBox.Text = avgTidalStrain.ToString("F4");

            // Compute average strain rate (du/dx)
            double strainRate = 0.0;
            for (int i = 1; i < gridX - 1; i++)
            {
                for (int j = 0; j < gridZ; j++)
                {
                    strainRate += Math.Abs((uVelocityProfile[i + 1, j] - uVelocityProfile[i - 1, j]) / (2 * dx));
                }
            }
            strainRate /= (gridX - 2) * gridZ;
            strainRateTextBox.Text = strainRate.ToString("F4");

            // Compute average Richardson number and density ratio
            double avgRi = 0.0;
            double avgDensityRatio = 0.0;
            int count = 0;
            for (int i = 0; i < gridX; i++)
            {
                for (int j = 1; j < gridZ - 1; j++)
                {
                    double dudz = (uVelocityProfile[i, j + 1] - uVelocityProfile[i, j - 1]) / (2 * dz);
                    double dTdz = (temperatureProfile[i, j + 1] - temperatureProfile[i, j - 1]) / (2 * dz);
                    double dsdz = (salinityProfile[i, j + 1] - salinityProfile[i, j - 1]) / (2 * dz);
                    double dCdz = (sedimentProfile[i, j + 1] - sedimentProfile[i, j - 1]) / (2 * dz);
                    double drhodz = -rho0 * (alpha * dTdz - beta * dsdz - gamma * dCdz); // Density gradient
                    double N2 = -g / rho0 * drhodz; // Buoyancy frequency squared
                    double shear2 = dudz * dudz;
                    if (shear2 > 1e-10)
                    {
                        avgRi += N2 / shear2;
                        count++;
                    }
                    double alphaDeltaT = alpha * dTdz;
                    double betaDeltaS = beta * dsdz;
                    if (Math.Abs(betaDeltaS) > 1e-10)
                    {
                        avgDensityRatio += Math.Abs(alphaDeltaT / betaDeltaS);
                    }
                }
            }
            avgRi = count > 0 ? avgRi / count : 0.0;
            avgDensityRatio = count > 0 ? avgDensityRatio / count : 0.0;
            richardsonTextBox.Text = avgRi.ToString("F4");
            densityRatioTextBox.Text = avgDensityRatio.ToString("F4");
            richardsonHistory[historyIndex] = avgRi;
            densityRatioHistory[historyIndex] = avgDensityRatio;
            sedimentHistory[historyIndex] = sedimentProfile[gridX / 2, gridZ / 2]; // Mid-point
            historyIndex = (historyIndex + 1) % richardsonHistory.Length;

            timeTextBox.Text = currentTime.ToString("F2");
        }

        private void VisualizationBox_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            int width = visualizationBox.Width;
            int height = visualizationBox.Height;
            g.ScaleTransform(zoomLevel, zoomLevel);

            // Draw background
            g.FillRectangle(Brushes.LightBlue, 0, 0, width, height);

            string mode = visualizationModeComboBox.SelectedItem?.ToString();
            if (mode == "Salinity")
            {
                // Draw salinity profile
                for (int i = 0; i < gridX - 1; i++)
                {
                    for (int j = 0; j < gridZ - 1; j++)
                    {
                        int x = (int)(i * width / (gridX - 1) / zoomLevel);
                        int y = (int)(j * height / (gridZ - 1) / zoomLevel);
                        int cellWidth = (int)(width / (gridX - 1) / zoomLevel) + 1;
                        int cellHeight = (int)(height / (gridZ - 1) / zoomLevel) + 1;
                        float salinity = (float)salinityProfile[i, j];
                        int colorValue = (int)(salinity / 35.0 * 255);
                        Color cellColor = Color.FromArgb(255 - colorValue, 255 - colorValue, 255);
                        using (SolidBrush brush = new SolidBrush(cellColor))
                        {
                            g.FillRectangle(brush, x, y, cellWidth, cellHeight);
                        }
                    }
                }
                // Draw color bar
                for (int y = 0; y < height / 2; y++)
                {
                    int colorValue = (int)(y * 255.0 / (height / 2));
                    using (SolidBrush brush = new SolidBrush(Color.FromArgb(255 - colorValue, 255 - colorValue, 255)))
                    {
                        g.FillRectangle(brush, (int)(width / zoomLevel) - 20, y, 10, 1);
                    }
                }
                g.DrawString("0 PSU", new Font("Consolas", 8), Brushes.Black, (int)(width / zoomLevel) - 30, height / 2 - 10);
                g.DrawString("35 PSU", new Font("Consolas", 8), Brushes.Black, (int)(width / zoomLevel) - 30, 0);
            }
            else if (mode == "Vorticity")
            {
                // Draw vorticity (dw/dx - du/dz)
                double maxVorticity = 0.0;
                for (int i = 1; i < gridX - 1; i++)
                {
                    for (int j = 1; j < gridZ - 1; j++)
                    {
                        double dwdx = (wVelocityProfile[i + 1, j] - wVelocityProfile[i - 1, j]) / (2 * estuaryLength / (gridX - 1));
                        double dudz = (uVelocityProfile[i, j + 1] - uVelocityProfile[i, j - 1]) / (2 * estuaryDepth / (gridZ - 1));
                        maxVorticity = Math.Max(maxVorticity, Math.Abs(dwdx - dudz));
                    }
                }
                for (int i = 1; i < gridX - 1; i++)
                {
                    for (int j = 1; j < gridZ - 1; j++)
                    {
                        int x = (int)(i * width / (gridX - 1) / zoomLevel);
                        int y = (int)(j * height / (gridZ - 1) / zoomLevel);
                        int cellWidth = (int)(width / (gridX - 1) / zoomLevel) + 1;
                        int cellHeight = (int)(height / (gridZ - 1) / zoomLevel) + 1;
                        double dwdx = (wVelocityProfile[i + 1, j] - wVelocityProfile[i - 1, j]) / (2 * estuaryLength / (gridX - 1));
                        double dudz = (uVelocityProfile[i, j + 1] - uVelocityProfile[i, j - 1]) / (2 * estuaryDepth / (gridZ - 1));
                        double vorticity = dwdx - dudz;
                        int colorValue = maxVorticity > 0 ? (int)(Math.Abs(vorticity) / maxVorticity * 255) : 0;
                        Color cellColor = Color.FromArgb(colorValue, 0, 255 - colorValue);
                        using (SolidBrush brush = new SolidBrush(cellColor))
                        {
                            g.FillRectangle(brush, x, y, cellWidth, cellHeight);
                        }
                    }
                }
                // Draw color bar for vorticity
                for (int y = 0; y < height / 2; y++)
                {
                    int colorValue = (int)(y * 255.0 / (height / 2));
                    using (SolidBrush brush = new SolidBrush(Color.FromArgb(colorValue, 0, 255 - colorValue)))
                    {
                        g.FillRectangle(brush, (int)(width / zoomLevel) - 20, y, 10, 1);
                    }
                }
                g.DrawString("0", new Font("Consolas", 8), Brushes.Black, (int)(width / zoomLevel) - 30, height / 2 - 10);
                g.DrawString($"{maxVorticity:F2}", new Font("Consolas", 8), Brushes.Black, (int)(width / zoomLevel) - 30, 0);
            }
            else if (mode == "Density")
            {
                // Draw density profile
                double minDensity = double.MaxValue, maxDensity = double.MinValue;
                for (int i = 0; i < gridX; i++)
                {
                    for (int j = 0; j < gridZ; j++)
                    {
                        double density = ComputeDensity(i, j);
                        minDensity = Math.Min(minDensity, density);
                        maxDensity = Math.Max(maxDensity, density);
                    }
                }
                for (int i = 0; i < gridX - 1; i++)
                {
                    for (int j = 0; j < gridZ - 1; j++)
                    {
                        int x = (int)(i * width / (gridX - 1) / zoomLevel);
                        int y = (int)(j * height / (gridZ - 1) / zoomLevel);
                        int cellWidth = (int)(width / (gridX - 1) / zoomLevel) + 1;
                        int cellHeight = (int)(height / (gridZ - 1) / zoomLevel) + 1;
                        double density = ComputeDensity(i, j);
                        int colorValue = maxDensity > minDensity ?
                            (int)((density - minDensity) / (maxDensity - minDensity) * 255) : 0;
                        Color cellColor = Color.FromArgb(255 - colorValue, 255, 255 - colorValue);
                        using (SolidBrush brush = new SolidBrush(cellColor))
                        {
                            g.FillRectangle(brush, x, y, cellWidth, cellHeight);
                        }
                    }
                }
                // Draw color bar for density
                for (int y = 0; y < height / 2; y++)
                {
                    int colorValue = (int)(y * 255.0 / (height / 2));
                    using (SolidBrush brush = new SolidBrush(Color.FromArgb(255 - colorValue, 255, 255 - colorValue)))
                    {
                        g.FillRectangle(brush, (int)(width / zoomLevel) - 20, y, 10, 1);
                    }
                }
                g.DrawString($"{minDensity:F2}", new Font("Consolas", 8), Brushes.Black, (int)(width / zoomLevel) - 30, height / 2 - 10);
                g.DrawString($"{maxDensity:F2}", new Font("Consolas", 8), Brushes.Black, (int)(width / zoomLevel) - 30, 0);
            }
            else if (mode == "Turbidity")
            {
                // Draw suspended sediment (turbidity) profile
                double maxSediment = 0.0;
                for (int i = 0; i < gridX; i++)
                {
                    for (int j = 0; j < gridZ; j++)
                    {
                        maxSediment = Math.Max(maxSediment, sedimentProfile[i, j]);
                    }
                }
                for (int i = 0; i < gridX - 1; i++)
                {
                    for (int j = 0; j < gridZ - 1; j++)
                    {
                        int x = (int)(i * width / (gridX - 1) / zoomLevel);
                        int y = (int)(j * height / (gridZ - 1) / zoomLevel);
                        int cellWidth = (int)(width / (gridX - 1) / zoomLevel) + 1;
                        int cellHeight = (int)(height / (gridZ - 1) / zoomLevel) + 1;
                        double sediment = sedimentProfile[i, j];
                        int colorValue = maxSediment > 0 ? (int)(sediment / maxSediment * 255) : 0;
                        Color cellColor = Color.FromArgb(255, 255 - colorValue, 255 - colorValue);
                        using (SolidBrush brush = new SolidBrush(cellColor))
                        {
                            g.FillRectangle(brush, x, y, cellWidth, cellHeight);
                        }
                    }
                }
                // Draw color bar for turbidity
                for (int y = 0; y < height / 2; y++)
                {
                    int colorValue = (int)(y * 255.0 / (height / 2));
                    using (SolidBrush brush = new SolidBrush(Color.FromArgb(255, 255 - colorValue, 255 - colorValue)))
                    {
                        g.FillRectangle(brush, (int)(width / zoomLevel) - 20, y, 10, 1);
                    }
                }
                g.DrawString("0 g/m³", new Font("Consolas", 8), Brushes.Black, (int)(width / zoomLevel) - 30, height / 2 - 10);
                g.DrawString($"{maxSediment:F2}", new Font("Consolas", 8), Brushes.Black, (int)(width / zoomLevel) - 30, 0);
            }
            else if (mode == "Bedload")
            {
                // Draw bedload profile (1D along bottom)
                double maxBedload = 0.0;
                for (int i = 0; i < gridX; i++)
                {
                    maxBedload = Math.Max(maxBedload, bedloadProfile[i]);
                }
                for (int i = 0; i < gridX - 1; i++)
                {
                    int x = (int)(i * width / (gridX - 1) / zoomLevel);
                    int cellWidth = (int)(width / (gridX - 1) / zoomLevel) + 1;
                    int cellHeight = (int)(height / zoomLevel / 10); // Thin strip at bottom
                    double bedload = bedloadProfile[i];
                    int colorValue = maxBedload > 0 ? (int)(bedload / maxBedload * 255) : 0;
                    Color cellColor = Color.FromArgb(139, 69, 19); // Brown for sediment
                    using (SolidBrush brush = new SolidBrush(Color.FromArgb(colorValue, cellColor)))
                    {
                        g.FillRectangle(brush, x, (int)(height / zoomLevel) - cellHeight, cellWidth, cellHeight);
                    }
                }
                // Draw color bar for bedload
                for (int y = 0; y < height / 2; y++)
                {
                    int colorValue = (int)(y * 255.0 / (height / 2));
                    using (SolidBrush brush = new SolidBrush(Color.FromArgb(colorValue, 139, 69, 19)))
                    {
                        g.FillRectangle(brush, (int)(width / zoomLevel) - 20, y, 10, 1);
                    }
                }
                g.DrawString("0 kg/m²", new Font("Consolas", 8), Brushes.Black, (int)(width / zoomLevel) - 30, height / 2 - 10);
                g.DrawString($"{maxBedload:F2}", new Font("Consolas", 8), Brushes.Black, (int)(width / zoomLevel) - 30, 0);
            }

            // Draw velocity vectors
            for (int i = 0; i < gridX; i += 20)
            {
                for (int j = 0; j < gridZ; j += 10)
                {
                    int x = (int)(i * width / (gridX - 1) / zoomLevel);
                    int y = (int)(j * height / (gridZ - 1) / zoomLevel);
                    int uLength = (int)(uVelocityProfile[i, j] * 50 / zoomLevel);
                    int wLength = (int)(wVelocityProfile[i, j] * 50 / zoomLevel);
                    g.DrawLine(Pens.Red, x, y, x + uLength, y + wLength);
                }
            }

            // Draw tidal phase indicator
            double tidalPhase = 2 * Math.PI * currentTime / tidalPeriod;
            string phaseText = (Math.Cos(tidalPhase) > 0) ? "Flood" : "Ebb";
            g.DrawString($"Tidal Phase: {phaseText}", new Font("Consolas", 10), Brushes.Black, 10 / zoomLevel, 10 / zoomLevel);
        }

        private void TimeSeriesBox_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            int width = timeSeriesBox.Width;
            int height = timeSeriesBox.Height;
            g.FillRectangle(Brushes.White, 0, 0, width, height);

            // Draw axes
            g.DrawLine(Pens.Black, 50, height - 30, width - 10, height - 30);
            g.DrawLine(Pens.Black, 50, 10, 50, height - 30);
            g.DrawString("Time", new Font("Consolas", 8), Brushes.Black, width - 40, height - 20);
            g.DrawString("Ri / Rρ / C", new Font("Consolas", 8), Brushes.Black, 30, 10);

            // Find max values for scaling
            double maxRi = 0.0, maxDensityRatio = 0.0, maxSediment = 0.0;
            for (int i = 0; i < richardsonHistory.Length; i++)
            {
                maxRi = Math.Max(maxRi, Math.Abs(richardsonHistory[i]));
                maxDensityRatio = Math.Max(maxDensityRatio, Math.Abs(densityRatioHistory[i]));
                maxSediment = Math.Max(maxSediment, sedimentHistory[i]);
            }
            maxRi = Math.Max(maxRi, 1.0);
            maxDensityRatio = Math.Max(maxDensityRatio, 1.0);
            maxSediment = Math.Max(maxSediment, 100.0);

            // Plot Ri history
            for (int i = 1; i < historyIndex; i++)
            {
                int x1 = 50 + (i - 1) * (width - 60) / richardsonHistory.Length;
                int x2 = 50 + i * (width - 60) / richardsonHistory.Length;
                int y1 = (int)(height - 30 - (richardsonHistory[i - 1] / maxRi) * (height - 40) / 3);
                int y2 = (int)(height - 30 - (richardsonHistory[i] / maxRi) * (height - 40) / 3);
                g.DrawLine(Pens.Blue, x1, y1, x2, y2);
            }

            // Plot density ratio history
            for (int i = 1; i < historyIndex; i++)
            {
                int x1 = 50 + (i - 1) * (width - 60) / densityRatioHistory.Length;
                int x2 = 50 + i * (width - 60) / densityRatioHistory.Length;
                int y1 = (int)(height - 30 - (densityRatioHistory[i - 1] / maxDensityRatio) * (height - 40) / 3);
                int y2 = (int)(height - 30 - (densityRatioHistory[i] / maxDensityRatio) * (height - 40) / 3);
                g.DrawLine(Pens.Red, x1, y1, x2, y2);
            }

            // Plot sediment concentration history
            for (int i = 1; i < historyIndex; i++)
            {
                int x1 = 50 + (i - 1) * (width - 60) / sedimentHistory.Length;
                int x2 = 50 + i * (width - 60) / sedimentHistory.Length;
                int y1 = (int)(height - 30 - (sedimentHistory[i - 1] / maxSediment) * (height - 40) / 3);
                int y2 = (int)(height - 30 - (sedimentHistory[i] / maxSediment) * (height - 40) / 3);
                g.DrawLine(Pens.Green, x1, y1, x2, y2);
            }

            // Draw legend
            g.DrawLine(Pens.Blue, 60, 20, 80, 20);
            g.DrawString("Ri", new Font("Consolas", 8), Brushes.Black, 90, 15);
            g.DrawLine(Pens.Red, 60, 35, 80, 35);
            g.DrawString("Rρ", new Font("Consolas", 8), Brushes.Black, 90, 30);
            g.DrawLine(Pens.Green, 60, 50, 80, 50);
            g.DrawString("C", new Font("Consolas", 8), Brushes.Black, 90, 45);
        }

        private void StartButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Start();
            isSimulationRunning = true;
            UpdateButtonStates();
        }

        private void PauseButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isSimulationRunning = false;
            UpdateButtonStates();
        }

        private void ResetButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isSimulationRunning = false;
            currentTime = 0.0;
            timeTextBox.Text = "0.0";
            strainRateTextBox.Text = "0.0";
            richardsonTextBox.Text = "0.0";
            tidalStrainTextBox.Text = "0.0";
            densityRatioTextBox.Text = "0.0";
            InitializeProfiles();
            Array.Clear(uVelocityProfile, 0, uVelocityProfile.Length);
            Array.Clear(wVelocityProfile, 0, wVelocityProfile.Length);
            Array.Clear(richardsonHistory, 0, richardsonHistory.Length);
            Array.Clear(densityRatioHistory, 0, densityRatioHistory.Length);
            Array.Clear(sedimentHistory, 0, sedimentHistory.Length);
            historyIndex = 0;
            visualizationBox.Invalidate();
            timeSeriesBox.Invalidate();
            UpdateButtonStates();
        }

        private void UpdateButtonStates()
        {
            startButton.Enabled = !isSimulationRunning;
            pauseButton.Enabled = isSimulationRunning;
            resetButton.Enabled = true;
        }

        public void ShowTidalStrainWindow()
        {
            if (tidalStrainForm == null || tidalStrainForm.IsDisposed)
            {
                InitializeForm(); // Recreate the form if it is null or disposed
            }
            if (!tidalStrainForm.Visible)
            {
                tidalStrainForm.Show();
            }
            else
            {
                tidalStrainForm.BringToFront();
            }
        }
    }
}