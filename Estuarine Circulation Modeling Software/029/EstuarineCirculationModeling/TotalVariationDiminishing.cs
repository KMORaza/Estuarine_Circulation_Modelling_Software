using System;
using System.Drawing;
using System.Drawing.Drawing2D;
using System.Windows.Forms;

namespace EstuarineCirculationModeling
{
    public class TotalVariationDiminishing
    {
        private Form tvdWindow;
        private Panel visualizationPanel;
        private TextBox outputTextBox;
        private Button startButton;
        private Button pauseButton;
        private Button resetButton;
        private RadioButton planViewRadio;
        private RadioButton crossSectionRadio;
        private RadioButton contourPlotRadio;
        private RadioButton quiverPlotRadio;
        private Label depthLabel;
        private TextBox depthTextBox;
        private Label xPositionLabel;
        private TextBox xPositionTextBox;
        private Label tidalPeriodLabel;
        private TextBox tidalPeriodTextBox;
        private Label turbulenceModelLabel;
        private ComboBox turbulenceModelComboBox;
        private Timer simulationTimer;
        private double[,] salinityField;
        private double[,] temperatureField;
        private double[,] velocityFieldX;
        private double[,] velocityFieldZ;
        private double[,] densityField;
        private double[,] kField;
        private double[,] epsilonField;
        private double[,] omegaField;
        private double[,] eddyViscosity;
        private double[,] eddyDiffusivity;
        private double dx, dz, dt, time;
        private int gridSizeX = 200;
        private int gridSizeZ = 100;
        private bool isRunning;
        private string visualizationMode = "PlanView";
        private string turbulenceModel = "k-epsilon";
        private double selectedDepth = 5.0;
        private double selectedXPosition = 500.0;
        private double estuaryLength = 1000.0;
        private double estuaryDepth = 10.0;
        private double tidalPeriod = 60.0;
        private double tidalAmplitude = 0.3;
        private double horizontalDiffusion = 0.01;
        private double g = 9.81;
        private double rho0 = 1000.0;
        private double alpha = 2e-4;
        private double beta_s = 8e-4;
        private double T0 = 20.0;
        private double S0 = 35.0;
        // k-ε model constants
        private double C_mu = 0.09;
        private double sigma_k = 1.0;
        private double sigma_epsilon = 1.3;
        private double C1_epsilon = 1.44;
        private double C2_epsilon = 1.92;
        // k-ω model constants
        private double sigma_k_omega = 0.5;
        private double sigma_omega = 0.5;
        private double beta_star = 0.09;
        private double beta = 0.075;
        private double alpha_omega = 0.52;
        private double riverSalinitySource = 0.0;
        private double oceanSalinitySource = 0.1;
        private double riverTempSource = 0.05;
        private double oceanTempSource = 0.02;

        public TotalVariationDiminishing()
        {
            dx = estuaryLength / gridSizeX;
            dz = estuaryDepth / gridSizeZ;
            dt = 0.025;
            time = 0.0;
            isRunning = false;
            salinityField = new double[gridSizeX, gridSizeZ];
            temperatureField = new double[gridSizeX, gridSizeZ];
            velocityFieldX = new double[gridSizeX, gridSizeZ];
            velocityFieldZ = new double[gridSizeX, gridSizeZ];
            densityField = new double[gridSizeX, gridSizeZ];
            kField = new double[gridSizeX, gridSizeZ];
            epsilonField = new double[gridSizeX, gridSizeZ];
            omegaField = new double[gridSizeX, gridSizeZ];
            eddyViscosity = new double[gridSizeX, gridSizeZ];
            eddyDiffusivity = new double[gridSizeX, gridSizeZ];
            InitializeFields();
        }

        private void InitializeFields()
        {
            // Initialize salinity, temperature, and density with Boussinesq approximation
            for (int i = 0; i < gridSizeX; i++)
            {
                double x = i * dx;
                double x0 = 500.0;
                for (int j = 0; j < gridSizeZ; j++)
                {
                    double z = j * dz;
                    salinityField[i, j] = 17.5 * (1.0 + Math.Tanh((x - x0) / 50.0)) * (1.0 - 0.5 * z / estuaryDepth);
                    temperatureField[i, j] = Math.Max(15.0, 20.0 + 2.0 * (1.0 - Math.Tanh((x - x0) / 50.0)) - 2.0 * z / estuaryDepth);
                    // Boussinesq: Density variations for buoyancy only
                    densityField[i, j] = rho0 * (1.0 - alpha * (temperatureField[i, j] - T0) + beta_s * (salinityField[i, j] - S0));
                    densityField[i, j] = Math.Max(950.0, Math.Min(1050.0, densityField[i, j])); // Clamp for stability
                    velocityFieldX[i, j] = tidalAmplitude * (1.0 - z / estuaryDepth);
                    velocityFieldZ[i, j] = -0.02 * (z / estuaryDepth);
                    kField[i, j] = 1e-4;
                    epsilonField[i, j] = C_mu * Math.Pow(kField[i, j], 1.5) / 0.1;
                    omegaField[i, j] = kField[i, j] / (0.09 * 0.1);
                    eddyViscosity[i, j] = C_mu * kField[i, j] * kField[i, j] / epsilonField[i, j];
                    eddyDiffusivity[i, j] = eddyViscosity[i, j] / 0.7;
                }
            }
        }

        private void UpdateVelocityField()
        {
            // Boussinesq approximation: Density affects only buoyancy term
            double tidalPhase = 2.0 * Math.PI * time / tidalPeriod;
            for (int i = 0; i < gridSizeX; i++)
            {
                int i_prev = (i == 0) ? gridSizeX - 1 : i - 1;
                int i_next = (i == gridSizeX - 1) ? 0 : i + 1;
                for (int j = 0; j < gridSizeZ; j++)
                {
                    double z = j * dz;
                    // Tidal velocity in x-direction
                    velocityFieldX[i, j] = tidalAmplitude * Math.Sin(tidalPhase) * (1.0 - z / estuaryDepth);
                    // Buoyancy term for vertical velocity (Boussinesq)
                    double drho_dx = (i == 0 || i == gridSizeX - 1) ? 0 : (densityField[i_next, j] - densityField[i_prev, j]) / (2 * dx);
                    double buoyancy = -g * drho_dx / rho0;
                    velocityFieldZ[i, j] += dt * buoyancy;
                    velocityFieldZ[i, j] = Math.Max(-0.5, Math.Min(0.5, velocityFieldZ[i, j])); // Clamp for stability
                }
            }
        }

        public void ShowTVDWindow()
        {
            tvdWindow = new Form
            {
                Text = "Total Variation Diminishing",
                Size = new Size(800, 750),
                FormBorderStyle = FormBorderStyle.FixedDialog,
                MaximizeBox = false,
                StartPosition = FormStartPosition.CenterScreen
            };

            visualizationPanel = new Panel
            {
                Location = new Point(10, 10),
                Size = new Size(760, 400),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            visualizationPanel.Paint += VisualizationPanel_Paint;

            outputTextBox = new TextBox
            {
                Location = new Point(10, 420),
                Size = new Size(760, 100),
                Multiline = true,
                ReadOnly = true,
                ScrollBars = ScrollBars.Vertical,
                Font = new Font("Courier New", 10F),
                BackColor = Color.Black,
                ForeColor = Color.White,
                BorderStyle = BorderStyle.FixedSingle
            };

            planViewRadio = new RadioButton
            {
                Location = new Point(10, 530),
                Size = new Size(100, 25),
                Text = "Plan View",
                Checked = true,
                Font = new Font("Verdana", 8F)
            };
            planViewRadio.CheckedChanged += (s, e) => { if (planViewRadio.Checked) visualizationMode = "PlanView"; visualizationPanel.Invalidate(); };

            crossSectionRadio = new RadioButton
            {
                Location = new Point(120, 530),
                Size = new Size(120, 25),
                Text = "Cross-Section",
                Font = new Font("Verdana", 8F)
            };
            crossSectionRadio.CheckedChanged += (s, e) => { if (crossSectionRadio.Checked) visualizationMode = "CrossSection"; visualizationPanel.Invalidate(); };

            contourPlotRadio = new RadioButton
            {
                Location = new Point(250, 530),
                Size = new Size(100, 25),
                Text = "Contour Plot",
                Font = new Font("Verdana", 8F)
            };
            contourPlotRadio.CheckedChanged += (s, e) => { if (contourPlotRadio.Checked) visualizationMode = "ContourPlot"; visualizationPanel.Invalidate(); };

            quiverPlotRadio = new RadioButton
            {
                Location = new Point(360, 530),
                Size = new Size(100, 25),
                Text = "Quiver Plot",
                Font = new Font("Verdana", 8F)
            };
            quiverPlotRadio.CheckedChanged += (s, e) => { if (quiverPlotRadio.Checked) visualizationMode = "QuiverPlot"; visualizationPanel.Invalidate(); };

            depthLabel = new Label
            {
                Location = new Point(10, 560),
                Size = new Size(80, 25),
                Text = "Depth (m):",
                Font = new Font("Verdana", 8F)
            };

            depthTextBox = new TextBox
            {
                Location = new Point(90, 560),
                Size = new Size(100, 25),
                Text = selectedDepth.ToString(),
                Font = new Font("Verdana", 8F)
            };
            depthTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(depthTextBox.Text, out double depth))
                {
                    selectedDepth = Math.Max(0, Math.Min(estuaryDepth, depth));
                    visualizationPanel.Invalidate();
                }
            };

            xPositionLabel = new Label
            {
                Location = new Point(200, 560),
                Size = new Size(80, 25),
                Text = "X-Position (m):",
                Font = new Font("Verdana", 8F)
            };

            xPositionTextBox = new TextBox
            {
                Location = new Point(280, 560),
                Size = new Size(100, 25),
                Text = selectedXPosition.ToString(),
                Font = new Font("Verdana", 8F)
            };
            xPositionTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(xPositionTextBox.Text, out double xPos))
                {
                    selectedXPosition = Math.Max(0, Math.Min(estuaryLength, xPos));
                    visualizationPanel.Invalidate();
                }
            };

            tidalPeriodLabel = new Label
            {
                Location = new Point(390, 560),
                Size = new Size(80, 25),
                Text = "Tidal Period (s):",
                Font = new Font("Verdana", 8F)
            };

            tidalPeriodTextBox = new TextBox
            {
                Location = new Point(470, 560),
                Size = new Size(100, 25),
                Text = tidalPeriod.ToString(),
                Font = new Font("Verdana", 8F)
            };
            tidalPeriodTextBox.TextChanged += (s, e) =>
            {
                if (double.TryParse(tidalPeriodTextBox.Text, out double period))
                {
                    tidalPeriod = Math.Max(10.0, period);
                    visualizationPanel.Invalidate();
                }
            };

            turbulenceModelLabel = new Label
            {
                Location = new Point(10, 590),
                Size = new Size(100, 25),
                Text = "Turbulence Model:",
                Font = new Font("Verdana", 8F)
            };

            turbulenceModelComboBox = new ComboBox
            {
                Location = new Point(110, 590),
                Size = new Size(100, 25),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = new Font("Verdana", 8F)
            };
            turbulenceModelComboBox.Items.AddRange(new object[] { "k-epsilon", "k-omega" });
            turbulenceModelComboBox.SelectedIndex = 0;
            turbulenceModelComboBox.SelectedIndexChanged += (s, e) =>
            {
                turbulenceModel = turbulenceModelComboBox.SelectedItem.ToString();
                InitializeFields();
                visualizationPanel.Invalidate();
            };

            startButton = new Button
            {
                Location = new Point(220, 590),
                Size = new Size(100, 25),
                Text = "Start",
                FlatStyle = FlatStyle.Flat,
                Font = new Font("Verdana", 8F)
            };
            startButton.Click += StartButton_Click;

            pauseButton = new Button
            {
                Location = new Point(330, 590),
                Size = new Size(100, 25),
                Text = "Pause",
                FlatStyle = FlatStyle.Flat,
                Font = new Font("Verdana", 8F),
                Enabled = false
            };
            pauseButton.Click += PauseButton_Click;

            resetButton = new Button
            {
                Location = new Point(440, 590),
                Size = new Size(100, 25),
                Text = "Reset",
                FlatStyle = FlatStyle.Flat,
                Font = new Font("Verdana", 8F)
            };
            resetButton.Click += ResetButton_Click;

            tvdWindow.Controls.Add(visualizationPanel);
            tvdWindow.Controls.Add(outputTextBox);
            tvdWindow.Controls.Add(planViewRadio);
            tvdWindow.Controls.Add(crossSectionRadio);
            tvdWindow.Controls.Add(contourPlotRadio);
            tvdWindow.Controls.Add(quiverPlotRadio);
            tvdWindow.Controls.Add(depthLabel);
            tvdWindow.Controls.Add(depthTextBox);
            tvdWindow.Controls.Add(xPositionLabel);
            tvdWindow.Controls.Add(xPositionTextBox);
            tvdWindow.Controls.Add(tidalPeriodLabel);
            tvdWindow.Controls.Add(tidalPeriodTextBox);
            tvdWindow.Controls.Add(turbulenceModelLabel);
            tvdWindow.Controls.Add(turbulenceModelComboBox);
            tvdWindow.Controls.Add(startButton);
            tvdWindow.Controls.Add(pauseButton);
            tvdWindow.Controls.Add(resetButton);

            simulationTimer = new Timer
            {
                Interval = 100
            };
            simulationTimer.Tick += (s, e) => UpdateSimulation();

            UpdateButtonStates();
            tvdWindow.Show();
        }

        private double HLLFlux(double uL, double uR, double v, double dx, double dt)
        {
            double sL = v - Math.Abs(v);
            double sR = v + Math.Abs(v);
            if (sL >= 0)
                return v * uL;
            else if (sR <= 0)
                return v * uR;
            else
                return (sR * v * uL - sL * v * uR + sL * sR * (uR - uL)) / (sR - sL);
        }

        private void UpdateSimulation()
        {
            UpdateVelocityField();
            double[,] newSalinity = new double[gridSizeX, gridSizeZ];
            double[,] newTemperature = new double[gridSizeX, gridSizeZ];
            double[,] newDensity = new double[gridSizeX, gridSizeZ];
            double[,] newK = new double[gridSizeX, gridSizeZ];
            double[,] newEpsilon = new double[gridSizeX, gridSizeZ];
            double[,] newOmega = new double[gridSizeX, gridSizeZ];
            double[,] newVelocityZ = new double[gridSizeX, gridSizeZ];

            for (int i = 0; i < gridSizeX; i++)
            {
                int i_prev = (i == 0) ? gridSizeX - 1 : i - 1;
                int i_next = (i == gridSizeX - 1) ? 0 : i + 1;
                for (int j = 0; j < gridSizeZ; j++)
                {
                    int j_prev = (j == 0) ? 0 : j - 1; // No-flux at surface
                    int j_next = (j == gridSizeZ - 1) ? gridSizeZ - 1 : j + 1; // No-flux at bottom
                    double x = i * dx;

                    // Source terms
                    double salinitySource = x < 200.0 ? riverSalinitySource : (x > 800.0 ? oceanSalinitySource : 0.0);
                    double tempSource = x < 200.0 ? riverTempSource : (x > 800.0 ? oceanTempSource : 0.0);

                    // Shear and buoyancy production (Boussinesq)
                    double dudz = (j == 0 ? (velocityFieldX[i, 1] - velocityFieldX[i, 0]) / dz :
                                  j == gridSizeZ - 1 ? (velocityFieldX[i, gridSizeZ - 1] - velocityFieldX[i, gridSizeZ - 2]) / dz :
                                  (velocityFieldX[i, j + 1] - velocityFieldX[i, j - 1]) / (2 * dz));
                    double drho_dz = (j == 0 || j == gridSizeZ - 1) ? 0 : (densityField[i, j_next] - densityField[i, j_prev]) / (2 * dz);
                    double P_shear = eddyViscosity[i, j] * dudz * dudz;
                    double P_buoy = -g * eddyDiffusivity[i, j] * drho_dz / rho0;
                    double P = Math.Max(P_shear + P_buoy, 0);

                    // Salinity advection (TVD-HLL)
                    double r_sx = (salinityField[i, j] - salinityField[i_prev, j]) / (salinityField[i_next, j] - salinityField[i, j] + 1e-10);
                    double flux_sx;
                    if (Math.Abs(r_sx) > 2.0 || Math.Abs(r_sx) < 0.5)
                    {
                        flux_sx = HLLFlux(salinityField[i_prev, j], salinityField[i, j], velocityFieldX[i, j], dx, dt);
                    }
                    else
                    {
                        double phi_sx = Math.Max(0, Math.Min(2 * r_sx, Math.Min(1, (1 + r_sx) / 2)));
                        double fluxLow_sx = velocityFieldX[i, j] * salinityField[i_prev, j];
                        double fluxHigh_sx = velocityFieldX[i, j] * (salinityField[i, j] - salinityField[i_prev, j]);
                        flux_sx = fluxLow_sx + phi_sx * (fluxHigh_sx - fluxLow_sx) / 2;
                    }
                    double xTerm_s = -dt / dx * (flux_sx - (i == 0 ? HLLFlux(salinityField[gridSizeX - 1, j], salinityField[0, j], velocityFieldX[gridSizeX - 1, j], dx, dt) : HLLFlux(salinityField[i_prev, j], salinityField[i, j], velocityFieldX[i_prev, j], dx, dt)));

                    double r_sz = (salinityField[i, j] - salinityField[i, j_prev]) / (salinityField[i, j_next] - salinityField[i, j] + 1e-10);
                    double flux_sz;
                    if (Math.Abs(r_sz) > 2.0 || Math.Abs(r_sz) < 0.5)
                    {
                        flux_sz = HLLFlux(salinityField[i, j_prev], salinityField[i, j], velocityFieldZ[i, j], dz, dt);
                    }
                    else
                    {
                        double phi_sz = Math.Max(0, Math.Min(2 * r_sz, Math.Min(1, (1 + r_sz) / 2)));
                        double fluxLow_sz = velocityFieldZ[i, j] * salinityField[i, j_prev];
                        double fluxHigh_sz = velocityFieldZ[i, j] * (salinityField[i, j] - salinityField[i, j_prev]);
                        flux_sz = fluxLow_sz + phi_sz * (fluxHigh_sz - fluxLow_sz) / 2;
                    }
                    double zTerm_s = -dt / dz * (flux_sz - (j == 0 ? 0 : HLLFlux(salinityField[i, j_prev], salinityField[i, j], velocityFieldZ[i, j_prev], dz, dt)));

                    double diff_s = horizontalDiffusion * (salinityField[i_next, j] - 2 * salinityField[i, j] + salinityField[i_prev, j]) / (dx * dx) +
                                    (j == 0 || j == gridSizeZ - 1 ? 0 : eddyDiffusivity[i, j] * (salinityField[i, j_next] - 2 * salinityField[i, j] + salinityField[i, j_prev]) / (dz * dz));

                    newSalinity[i, j] = Math.Max(0.0, salinityField[i, j] + xTerm_s + zTerm_s + dt * diff_s + dt * salinitySource);

                    // Temperature advection (TVD-HLL)
                    double r_tx = (temperatureField[i, j] - temperatureField[i_prev, j]) / (temperatureField[i_next, j] - temperatureField[i, j] + 1e-10);
                    double flux_tx;
                    if (Math.Abs(r_tx) > 2.0 || Math.Abs(r_tx) < 0.5)
                    {
                        flux_tx = HLLFlux(temperatureField[i_prev, j], temperatureField[i, j], velocityFieldX[i, j], dx, dt);
                    }
                    else
                    {
                        double phi_tx = Math.Max(0, Math.Min(2 * r_tx, Math.Min(1, (1 + r_tx) / 2)));
                        double fluxLow_tx = velocityFieldX[i, j] * temperatureField[i_prev, j];
                        double fluxHigh_tx = velocityFieldX[i, j] * (temperatureField[i, j] - temperatureField[i_prev, j]);
                        flux_tx = fluxLow_tx + phi_tx * (fluxHigh_tx - fluxLow_tx) / 2;
                    }
                    double xTerm_t = -dt / dx * (flux_tx - (i == 0 ? HLLFlux(temperatureField[gridSizeX - 1, j], temperatureField[0, j], velocityFieldX[gridSizeX - 1, j], dx, dt) : HLLFlux(temperatureField[i_prev, j], temperatureField[i, j], velocityFieldX[i_prev, j], dx, dt)));

                    double r_tz = (temperatureField[i, j] - temperatureField[i, j_prev]) / (temperatureField[i, j_next] - temperatureField[i, j] + 1e-10);
                    double flux_tz;
                    if (Math.Abs(r_tz) > 2.0 || Math.Abs(r_tz) < 0.5)
                    {
                        flux_tz = HLLFlux(temperatureField[i, j_prev], temperatureField[i, j], velocityFieldZ[i, j], dz, dt);
                    }
                    else
                    {
                        double phi_tz = Math.Max(0, Math.Min(2 * r_tz, Math.Min(1, (1 + r_tz) / 2)));
                        double fluxLow_tz = velocityFieldZ[i, j] * temperatureField[i, j_prev];
                        double fluxHigh_tz = velocityFieldZ[i, j] * (temperatureField[i, j] - temperatureField[i, j_prev]);
                        flux_tz = fluxLow_tz + phi_tz * (fluxHigh_tz - fluxLow_tz) / 2;
                    }
                    double zTerm_t = -dt / dz * (flux_tz - (j == 0 ? 0 : HLLFlux(temperatureField[i, j_prev], temperatureField[i, j], velocityFieldZ[i, j_prev], dz, dt)));

                    double diff_t = horizontalDiffusion * (temperatureField[i_next, j] - 2 * temperatureField[i, j] + temperatureField[i_prev, j]) / (dx * dx) +
                                    (j == 0 || j == gridSizeZ - 1 ? 0 : eddyDiffusivity[i, j] * (temperatureField[i, j_next] - 2 * temperatureField[i, j] + temperatureField[i, j_prev]) / (dz * dz));

                    newTemperature[i, j] = Math.Max(15.0, temperatureField[i, j] + xTerm_t + zTerm_t + dt * diff_t + dt * tempSource);

                    // Boussinesq: Update density for buoyancy
                    newDensity[i, j] = rho0 * (1.0 - alpha * (newTemperature[i, j] - T0) + beta_s * (newSalinity[i, j] - S0));
                    newDensity[i, j] = Math.Max(950.0, Math.Min(1050.0, newDensity[i, j])); // Clamp for Boussinesq stability

                    // Turbulence: k equation
                    double r_kx = (kField[i, j] - kField[i_prev, j]) / (kField[i_next, j] - kField[i, j] + 1e-10);
                    double flux_kx;
                    if (Math.Abs(r_kx) > 2.0 || Math.Abs(r_kx) < 0.5)
                    {
                        flux_kx = HLLFlux(kField[i_prev, j], kField[i, j], velocityFieldX[i, j], dx, dt);
                    }
                    else
                    {
                        double phi_kx = Math.Max(0, Math.Min(2 * r_kx, Math.Min(1, (1 + r_kx) / 2)));
                        double fluxLow_kx = velocityFieldX[i, j] * kField[i_prev, j];
                        double fluxHigh_kx = velocityFieldX[i, j] * (kField[i, j] - kField[i_prev, j]);
                        flux_kx = fluxLow_kx + phi_kx * (fluxHigh_kx - fluxLow_kx) / 2;
                    }
                    double xTerm_k = -dt / dx * (flux_kx - (i == 0 ? HLLFlux(kField[gridSizeX - 1, j], kField[0, j], velocityFieldX[gridSizeX - 1, j], dx, dt) : HLLFlux(kField[i_prev, j], kField[i, j], velocityFieldX[i_prev, j], dx, dt)));

                    double r_kz = (kField[i, j] - kField[i, j_prev]) / (kField[i, j_next] - kField[i, j] + 1e-10);
                    double flux_kz;
                    if (Math.Abs(r_kz) > 2.0 || Math.Abs(r_kz) < 0.5)
                    {
                        flux_kz = HLLFlux(kField[i, j_prev], kField[i, j], velocityFieldZ[i, j], dz, dt);
                    }
                    else
                    {
                        double phi_kz = Math.Max(0, Math.Min(2 * r_kz, Math.Min(1, (1 + r_kz) / 2)));
                        double fluxLow_kz = velocityFieldZ[i, j] * kField[i, j_prev];
                        double fluxHigh_kz = velocityFieldZ[i, j] * (kField[i, j] - kField[i, j_prev]);
                        flux_kz = fluxLow_kz + phi_kz * (fluxHigh_kz - fluxLow_kz) / 2;
                    }
                    double zTerm_k = -dt / dz * (flux_kz - (j == 0 ? 0 : HLLFlux(kField[i, j_prev], kField[i, j], velocityFieldZ[i, j_prev], dz, dt)));

                    double diff_k = (j == 0 || j == gridSizeZ - 1 ? 0 : (eddyViscosity[i, j] / sigma_k) * (kField[i, j_next] - 2 * kField[i, j] + kField[i, j_prev]) / (dz * dz));
                    double kSource = turbulenceModel == "k-epsilon" ? P - epsilonField[i, j] : P - beta_star * kField[i, j] * omegaField[i, j];
                    newK[i, j] = Math.Max(1e-6, kField[i, j] + xTerm_k + zTerm_k + dt * diff_k + dt * kSource);

                    // Turbulence: ε or ω equation
                    if (turbulenceModel == "k-epsilon")
                    {
                        double r_ex = (epsilonField[i, j] - epsilonField[i_prev, j]) / (epsilonField[i_next, j] - epsilonField[i, j] + 1e-10);
                        double flux_ex;
                        if (Math.Abs(r_ex) > 2.0 || Math.Abs(r_ex) < 0.5)
                        {
                            flux_ex = HLLFlux(epsilonField[i_prev, j], epsilonField[i, j], velocityFieldX[i, j], dx, dt);
                        }
                        else
                        {
                            double phi_ex = Math.Max(0, Math.Min(2 * r_ex, Math.Min(1, (1 + r_ex) / 2)));
                            double fluxLow_ex = velocityFieldX[i, j] * epsilonField[i_prev, j];
                            double fluxHigh_ex = velocityFieldX[i, j] * (epsilonField[i, j] - epsilonField[i_prev, j]);
                            flux_ex = fluxLow_ex + phi_ex * (fluxHigh_ex - fluxLow_ex) / 2;
                        }
                        double xTerm_e = -dt / dx * (flux_ex - (i == 0 ? HLLFlux(epsilonField[gridSizeX - 1, j], epsilonField[0, j], velocityFieldX[gridSizeX - 1, j], dx, dt) : HLLFlux(epsilonField[i_prev, j], epsilonField[i, j], velocityFieldX[i_prev, j], dx, dt)));

                        double r_ez = (epsilonField[i, j] - epsilonField[i, j_prev]) / (epsilonField[i, j_next] - epsilonField[i, j] + 1e-10);
                        double flux_ez;
                        if (Math.Abs(r_ez) > 2.0 || Math.Abs(r_ez) < 0.5)
                        {
                            flux_ez = HLLFlux(epsilonField[i, j_prev], epsilonField[i, j], velocityFieldZ[i, j], dz, dt);
                        }
                        else
                        {
                            double phi_ez = Math.Max(0, Math.Min(2 * r_ez, Math.Min(1, (1 + r_ez) / 2)));
                            double fluxLow_ez = velocityFieldZ[i, j] * epsilonField[i, j_prev];
                            double fluxHigh_ez = velocityFieldZ[i, j] * (epsilonField[i, j] - epsilonField[i, j_prev]);
                            flux_ez = fluxLow_ez + phi_ez * (fluxHigh_ez - fluxLow_ez) / 2;
                        }
                        double zTerm_e = -dt / dz * (flux_ez - (j == 0 ? 0 : HLLFlux(epsilonField[i, j_prev], epsilonField[i, j], velocityFieldZ[i, j_prev], dz, dt)));

                        double diff_e = (j == 0 || j == gridSizeZ - 1 ? 0 : (eddyViscosity[i, j] / sigma_epsilon) * (epsilonField[i, j_next] - 2 * epsilonField[i, j] + epsilonField[i, j_prev]) / (dz * dz));
                        double epsilonSource = C1_epsilon * (P * epsilonField[i, j] / kField[i, j]) - C2_epsilon * (epsilonField[i, j] * epsilonField[i, j] / kField[i, j]);
                        newEpsilon[i, j] = Math.Max(1e-6, epsilonField[i, j] + xTerm_e + zTerm_e + dt * diff_e + dt * epsilonSource);
                    }
                    else
                    {
                        double r_ox = (omegaField[i, j] - omegaField[i_prev, j]) / (omegaField[i_next, j] - omegaField[i, j] + 1e-10);
                        double flux_ox;
                        if (Math.Abs(r_ox) > 2.0 || Math.Abs(r_ox) < 0.5)
                        {
                            flux_ox = HLLFlux(omegaField[i_prev, j], omegaField[i, j], velocityFieldX[i, j], dx, dt);
                        }
                        else
                        {
                            double phi_ox = Math.Max(0, Math.Min(2 * r_ox, Math.Min(1, (1 + r_ox) / 2)));
                            double fluxLow_ox = velocityFieldX[i, j] * omegaField[i_prev, j];
                            double fluxHigh_ox = velocityFieldX[i, j] * (omegaField[i, j] - omegaField[i_prev, j]);
                            flux_ox = fluxLow_ox + phi_ox * (fluxHigh_ox - fluxLow_ox) / 2;
                        }
                        double xTerm_o = -dt / dx * (flux_ox - (i == 0 ? HLLFlux(omegaField[gridSizeX - 1, j], omegaField[0, j], velocityFieldX[gridSizeX - 1, j], dx, dt) : HLLFlux(omegaField[i_prev, j], omegaField[i, j], velocityFieldX[i_prev, j], dx, dt)));

                        double r_oz = (omegaField[i, j] - omegaField[i, j_prev]) / (omegaField[i, j_next] - omegaField[i, j] + 1e-10);
                        double flux_oz;
                        if (Math.Abs(r_oz) > 2.0 || Math.Abs(r_oz) < 0.5)
                        {
                            flux_oz = HLLFlux(omegaField[i, j_prev], omegaField[i, j], velocityFieldZ[i, j], dz, dt);
                        }
                        else
                        {
                            double phi_oz = Math.Max(0, Math.Min(2 * r_oz, Math.Min(1, (1 + r_oz) / 2)));
                            double fluxLow_oz = velocityFieldZ[i, j] * omegaField[i, j_prev];
                            double fluxHigh_oz = velocityFieldZ[i, j] * (omegaField[i, j] - omegaField[i, j_prev]);
                            flux_oz = fluxLow_oz + phi_oz * (fluxHigh_oz - fluxLow_oz) / 2;
                        }
                        double zTerm_o = -dt / dz * (flux_oz - (j == 0 ? 0 : HLLFlux(omegaField[i, j_prev], omegaField[i, j], velocityFieldZ[i, j_prev], dz, dt)));

                        double diff_o = (j == 0 || j == gridSizeZ - 1 ? 0 : (eddyViscosity[i, j] / sigma_omega) * (omegaField[i, j_next] - 2 * omegaField[i, j] + omegaField[i, j_prev]) / (dz * dz));
                        double omegaSource = alpha_omega * (P * omegaField[i, j] / kField[i, j]) - beta * omegaField[i, j] * omegaField[i, j];
                        newOmega[i, j] = Math.Max(1e-6, omegaField[i, j] + xTerm_o + zTerm_o + dt * diff_o + dt * omegaSource);
                    }

                    // Update eddy viscosity and diffusivity
                    eddyViscosity[i, j] = turbulenceModel == "k-epsilon" ?
                        C_mu * kField[i, j] * kField[i, j] / epsilonField[i, j] :
                        kField[i, j] / omegaField[i, j];
                    eddyDiffusivity[i, j] = eddyViscosity[i, j] / 0.7;

                    // Update vertical velocity with eddy viscosity
                    double diff_vz = (j == 0 || j == gridSizeZ - 1 ? 0 : eddyViscosity[i, j] * (velocityFieldZ[i, j_next] - 2 * velocityFieldZ[i, j] + velocityFieldZ[i, j_prev]) / (dz * dz));
                    newVelocityZ[i, j] = velocityFieldZ[i, j] + dt * diff_vz;
                    newVelocityZ[i, j] = Math.Max(-0.5, Math.Min(0.5, newVelocityZ[i, j])); // Clamp for stability
                }
            }

            salinityField = newSalinity;
            temperatureField = newTemperature;
            densityField = newDensity;
            kField = newK;
            if (turbulenceModel == "k-epsilon") epsilonField = newEpsilon;
            else omegaField = newOmega;
            velocityFieldZ = newVelocityZ;
            time += dt;
            visualizationPanel.Invalidate();
            UpdateOutputTextBox();
        }

        private void VisualizationPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            float panelWidth = visualizationPanel.Width;
            float panelHeight = visualizationPanel.Height;

            if (visualizationMode == "PlanView")
            {
                int j = (int)(selectedDepth / dz);
                j = Math.Max(0, Math.Min(gridSizeZ - 1, j));

                for (int i = 0; i < gridSizeX - 1; i++)
                {
                    float x1 = i * panelWidth / gridSizeX;
                    float x2 = (i + 1) * panelWidth / gridSizeX;
                    float y1 = panelHeight / 3 - (float)(salinityField[i, j] / 35.0 * panelHeight / 3);
                    float y2 = panelHeight / 3 - (float)(salinityField[i + 1, j] / 35.0 * panelHeight / 3);
                    g.DrawLine(new Pen(Color.Blue, 2), x1, y1, x2, y2);
                }

                for (int i = 0; i < gridSizeX - 1; i++)
                {
                    float x1 = i * panelWidth / gridSizeX;
                    float x2 = (i + 1) * panelWidth / gridSizeX;
                    float y1 = 2 * panelHeight / 3 - (float)((temperatureField[i, j] - 15.0) / 10.0 * panelHeight / 3);
                    float y2 = 2 * panelHeight / 3 - (float)((temperatureField[i + 1, j] - 15.0) / 10.0 * panelHeight / 3);
                    g.DrawLine(new Pen(Color.Red, 2), x1, y1, x2, y2);
                }

                for (int i = 0; i < gridSizeX - 1; i++)
                {
                    float x1 = i * panelWidth / gridSizeX;
                    float x2 = (i + 1) * panelWidth / gridSizeX;
                    double mag1 = Math.Sqrt(velocityFieldX[i, j] * velocityFieldX[i, j] + velocityFieldZ[i, j] * velocityFieldZ[i, j]);
                    double mag2 = Math.Sqrt(velocityFieldX[i + 1, j] * velocityFieldX[i + 1, j] + velocityFieldZ[i + 1, j] * velocityFieldZ[i + 1, j]);
                    float y1 = panelHeight - (float)(mag1 / tidalAmplitude * panelHeight / 3);
                    float y2 = panelHeight - (float)(mag2 / tidalAmplitude * panelHeight / 3);
                    g.DrawLine(new Pen(Color.Green, 2), x1, y1, x2, y2);
                }
            }
            else if (visualizationMode == "CrossSection")
            {
                int i = (int)(selectedXPosition / dx);
                i = Math.Max(0, Math.Min(gridSizeX - 1, i));

                for (int j = 0; j < gridSizeZ - 1; j++)
                {
                    float z1 = j * panelHeight / gridSizeZ;
                    float z2 = (j + 1) * panelHeight / gridSizeZ;
                    float x1 = panelWidth / 3 * (float)(salinityField[i, j] / 35.0);
                    float x2 = panelWidth / 3 * (float)(salinityField[i, j + 1] / 35.0);
                    g.DrawLine(new Pen(Color.Blue, 2), x1, z1, x2, z2);
                }

                for (int j = 0; j < gridSizeZ - 1; j++)
                {
                    float z1 = j * panelHeight / gridSizeZ;
                    float z2 = (j + 1) * panelHeight / gridSizeZ;
                    float x1 = panelWidth / 3 + panelWidth / 3 * (float)((temperatureField[i, j] - 15.0) / 10.0);
                    float x2 = panelWidth / 3 + panelWidth / 3 * (float)((temperatureField[i, j + 1] - 15.0) / 10.0);
                    g.DrawLine(new Pen(Color.Red, 2), x1, z1, x2, z2);
                }

                for (int j = 0; j < gridSizeZ - 1; j++)
                {
                    float z1 = j * panelHeight / gridSizeZ;
                    float z2 = (j + 1) * panelHeight / gridSizeZ;
                    double mag1 = Math.Sqrt(velocityFieldX[i, j] * velocityFieldX[i, j] + velocityFieldZ[i, j] * velocityFieldZ[i, j]);
                    double mag2 = Math.Sqrt(velocityFieldX[i, j + 1] * velocityFieldX[i, j + 1] + velocityFieldZ[i, j + 1] * velocityFieldZ[i, j + 1]);
                    float x1 = 2 * panelWidth / 3 + panelWidth / 3 * (float)(mag1 / tidalAmplitude);
                    float x2 = 2 * panelWidth / 3 + panelWidth / 3 * (float)(mag2 / tidalAmplitude);
                    g.DrawLine(new Pen(Color.Green, 2), x1, z1, x2, z2);
                }
            }
            else if (visualizationMode == "ContourPlot")
            {
                for (int i = 0; i < gridSizeX; i++)
                {
                    for (int j = 0; j < gridSizeZ; j++)
                    {
                        float x = i * panelWidth / gridSizeX;
                        float y = j * panelHeight / gridSizeZ;
                        float cellWidth = panelWidth / gridSizeX;
                        float cellHeight = panelHeight / gridSizeZ;

                        int salinityColor = (int)(salinityField[i, j] / 35.0 * 255);
                        salinityColor = Math.Max(0, Math.Min(255, salinityColor));
                        using (SolidBrush brush = new SolidBrush(Color.FromArgb(salinityColor, salinityColor, 255)))
                        {
                            g.FillRectangle(brush, x, y, cellWidth / 2, cellHeight);
                        }

                        int tempColor = (int)((temperatureField[i, j] - 15.0) / 10.0 * 255);
                        tempColor = Math.Max(0, Math.Min(255, tempColor));
                        using (SolidBrush brush = new SolidBrush(Color.FromArgb(255, tempColor, tempColor)))
                        {
                            g.FillRectangle(brush, x + panelWidth / 2, y, cellWidth / 2, cellHeight);
                        }
                    }
                }
            }
            else if (visualizationMode == "QuiverPlot")
            {
                int step = 5; // Plot every 5th point to avoid clutter
                float arrowScale = 1000.0f; // Scale factor for arrow length
                using (Pen arrowPen = new Pen(Color.Green, 1))
                {
                    arrowPen.CustomEndCap = new AdjustableArrowCap(3, 3);
                    for (int i = 0; i < gridSizeX; i += step)
                    {
                        for (int j = 0; j < gridSizeZ; j += step)
                        {
                            float x = i * panelWidth / gridSizeX;
                            float y = j * panelHeight / gridSizeZ;
                            float u = (float)(velocityFieldX[i, j] * arrowScale);
                            float v = (float)(-velocityFieldZ[i, j] * arrowScale); // Negative due to GDI+ y-axis
                            g.DrawLine(arrowPen, x, y, x + u, y + v);
                        }
                    }
                }
            }
        }

        private void UpdateOutputTextBox()
        {
            double maxSalinity = salinityField[0, 0];
            double maxTemperature = temperatureField[0, 0];
            double maxVelocity = Math.Sqrt(velocityFieldX[0, 0] * velocityFieldX[0, 0] + velocityFieldZ[0, 0] * velocityFieldZ[0, 0]);
            for (int i = 0; i < gridSizeX; i++)
            {
                for (int j = 0; j < gridSizeZ; j++)
                {
                    maxSalinity = Math.Max(maxSalinity, salinityField[i, j]);
                    maxTemperature = Math.Max(maxTemperature, temperatureField[i, j]);
                    double mag = Math.Sqrt(velocityFieldX[i, j] * velocityFieldX[i, j] + velocityFieldZ[i, j] * velocityFieldZ[i, j]);
                    maxVelocity = Math.Max(maxVelocity, mag);
                }
            }
            outputTextBox.AppendText($"Time: {time:F2}s | Max Salinity: {maxSalinity:F2} PSU | Max Temperature: {maxTemperature:F2} °C | Max Velocity: {maxVelocity:F4} m/s\r\n");
        }

        private void StartButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Start();
            isRunning = true;
            UpdateButtonStates();
            outputTextBox.AppendText("TVD simulation started.\r\n");
        }

        private void PauseButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isRunning = false;
            UpdateButtonStates();
            outputTextBox.AppendText("TVD simulation paused.\r\n");
        }

        private void ResetButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isRunning = false;
            time = 0.0;
            InitializeFields();
            visualizationPanel.Invalidate();
            outputTextBox.Clear();
            outputTextBox.AppendText("TVD simulation reset.\r\n");
            UpdateButtonStates();
        }

        private void UpdateButtonStates()
        {
            startButton.Enabled = !isRunning;
            pauseButton.Enabled = isRunning;
            resetButton.Enabled = true;
        }
    }
}