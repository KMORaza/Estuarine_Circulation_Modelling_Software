using System;
using System.Windows.Forms;
using System.Drawing;
using System.Linq;

namespace EstuarineCirculationModeling
{
    public class CompForcingMechanism : Form
    {
        private Panel controlPanel;
        private Panel visualizationPanel;
        private TextBox outputConsoleTextBox;
        private Button startButton, pauseButton, resetButton;
        private Label riverDischargeLabel, tidalAmplitudeLabel, tidalPeriodLabel, windSpeedLabel, windDirectionLabel, salinityGradientLabel, waveHeightLabel, stormSurgeAmplitudeLabel, solverTypeLabel, minDepthLabel;
        private TextBox riverDischargeTextBox, tidalAmplitudeTextBox, tidalPeriodTextBox, windSpeedTextBox, windDirectionTextBox, salinityGradientTextBox, waveHeightTextBox, stormSurgeAmplitudeTextBox, minDepthTextBox;
        private ComboBox solverTypeComboBox;
        private CheckBox enableWetDryCheckBox;
        private Timer simulationTimer;
        private bool isSimulationRunning;
        private string solverType;
        private ShallowWaterEq2D swe2DSolver;
        private WetAndDryAlgo wetDryAlgo;
        // 1D solver fields
        private double[] velocityProfile, salinityProfile, waterLevel;
        // 3D solver fields
        private double[,,] u, v, w, p, salinity, temperature;
        private double estuaryLength = 10000.0; // m
        private double estuaryWidth = 2000.0; // m
        private double estuaryDepth = 10.0; // m
        private double currentTime = 0.0; // s
        private double dt = 1.0; // Initial time step in seconds
        private int nx = 50, ny = 20, nz = 10; // 3D grid points
        private int gridPoints = 100; // 1D grid points
        private double dx, dy, dz; // 3D grid spacings
        private double riverDischarge = 0.1; // m³/s
        private double tidalAmplitude = 1.0; // m
        private double tidalPeriod = 43200.0; // s (12 hours)
        private double windSpeed = 5.0; // m/s
        private double windDirection = 0.0; // degrees
        private double salinityGradient = 0.0035; // PSU/m
        private double waveHeight = 0.5; // m
        private double stormSurgeAmplitude = 0.0; // m
        private double minDepth = 0.01; // m (default minimum depth for wetting/drying)
        private bool enableWetDry = false;
        private double g = 9.81; // m/s²
        private double rho0 = 1000.0; // kg/m³ (freshwater density)
        private double rhoOcean = 1025.0; // kg/m³ (ocean water density)
        private double nu = 1e-6; // Kinematic viscosity (m²/s)
        private double kappa = 1e-4; // Salinity diffusion coefficient (m²/s)
        private const double courantNumber = 0.4; // Stricter CFL safety factor
        private const double frictionCoefficient = 0.0025; // Bottom friction coefficient
        private const double waterLevelDamping = 1e-6; // Damping coefficient
        private const double coriolisParameter = 1e-4; // s^-1 (mid-latitude approximation)
        private const double eddyViscosity = 0.01; // m²/s (turbulence closure)
        private const double atmosphericPressureGradient = 0.0001; // Pa/m
        private const double seasonalSalinityAmplitude = 2.0; // PSU
        private WindForcing windForcing;
        private WaveCurrentInteraction waveInteraction; // Add WaveCurrentInteraction instance
        private double wavePeriod = 10.0; // Default wave period (s), adjust via UI if available

        public CompForcingMechanism()
        {
            InitializeComponents();
            InitializeSimulation();
            isSimulationRunning = false;
            UpdateButtonStates();
            this.FormClosing += CompForcingMechanism_FormClosing;
        }

        private void CompForcingMechanism_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (simulationTimer != null)
            {
                simulationTimer.Stop();
                simulationTimer.Dispose();
                isSimulationRunning = false;
            }
        }

        public static void ShowCompForcingWindow()
        {
            using (var form = new CompForcingMechanism())
            {
                form.ShowDialog();
            }
        }

        private void InitializeComponents()
        {
            this.Text = "Comprehensive Forcing";
            this.ClientSize = new Size(900, 700);
            this.FormBorderStyle = FormBorderStyle.FixedDialog;
            this.MaximizeBox = false;
            this.Font = new Font("Verdana", 8.25F);
            this.BackColor = SystemColors.Control;
            this.ForeColor = Color.Black;

            controlPanel = new Panel
            {
                Location = new Point(12, 12),
                Size = new Size(300, 520),
                TabIndex = 0,
                AutoScroll = true,
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = SystemColors.Control,
                ForeColor = Color.Black
            };

            int yOffset = 20;
            riverDischargeLabel = CreateLabel("River Discharge (m³/s):", 10, yOffset);
            riverDischargeTextBox = CreateTextBox("0.1", 10, yOffset + 20);
            yOffset += 50;
            tidalAmplitudeLabel = CreateLabel("Tidal Amplitude (m):", 10, yOffset);
            tidalAmplitudeTextBox = CreateTextBox("1.0", 10, yOffset + 20);
            yOffset += 50;
            tidalPeriodLabel = CreateLabel("Tidal Period (s):", 10, yOffset);
            tidalPeriodTextBox = CreateTextBox("43200", 10, yOffset + 20);
            yOffset += 50;
            windSpeedLabel = CreateLabel("Wind Speed (m/s):", 10, yOffset);
            windSpeedTextBox = CreateTextBox("5.0", 10, yOffset + 20);
            yOffset += 50;
            windDirectionLabel = CreateLabel("Wind Direction (deg):", 10, yOffset);
            windDirectionTextBox = CreateTextBox("0.0", 10, yOffset + 20);
            yOffset += 50;
            salinityGradientLabel = CreateLabel("Sal-Ein Grad (PSU/m):", 10, yOffset);
            salinityGradientTextBox = CreateTextBox("0.0035", 10, yOffset + 20);
            yOffset += 50;
            waveHeightLabel = CreateLabel("Wave Height (m):", 10, yOffset);
            waveHeightTextBox = CreateTextBox("0.5", 10, yOffset + 20);
            yOffset += 50;
            stormSurgeAmplitudeLabel = CreateLabel("Storm Surge Amp (m):", 10, yOffset);
            stormSurgeAmplitudeTextBox = CreateTextBox("0.0", 10, yOffset + 20);
            yOffset += 50;
            solverTypeLabel = CreateLabel("Solver Type:", 10, yOffset);
            solverTypeComboBox = new ComboBox
            {
                Location = new Point(10, yOffset + 20),
                Size = new Size(250, 22),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = new Font("Verdana", 8.25F),
                BackColor = Color.White,
                ForeColor = Color.Black
            };
            solverTypeComboBox.Items.AddRange(new[] { "1D Simplified", "2D Shallow Water", "3D Navier-Stokes" });
            solverTypeComboBox.SelectedIndex = 0; // Default to 1D
            solverTypeComboBox.SelectedIndexChanged += (s, e) => { solverType = solverTypeComboBox.SelectedItem.ToString(); InitializeSimulation(); };
            yOffset += 50;
            enableWetDryCheckBox = new CheckBox
            {
                Location = new Point(10, yOffset),
                Size = new Size(250, 22),
                Text = "Enable Wetting and Drying",
                Font = new Font("Verdana", 8.25F),
                BackColor = SystemColors.Control,
                ForeColor = Color.Black,
                Checked = false
            };
            enableWetDryCheckBox.CheckedChanged += (s, e) => { enableWetDry = enableWetDryCheckBox.Checked; InitializeSimulation(); };
            yOffset += 30;
            minDepthLabel = CreateLabel("Min Depth (m):", 10, yOffset);
            minDepthTextBox = CreateTextBox("0.01", 10, yOffset + 20);
            yOffset += 50;

            startButton = CreateButton("Start", 10, yOffset, startButton_Click);
            yOffset += 35;
            pauseButton = CreateButton("Pause", 10, yOffset, pauseButton_Click);
            yOffset += 35;
            resetButton = CreateButton("Reset", 10, yOffset, resetButton_Click);

            outputConsoleTextBox = new TextBox
            {
                Location = new Point(12, 540),
                Size = new Size(860, 150),
                Multiline = true,
                ReadOnly = true,
                ScrollBars = ScrollBars.Vertical,
                Font = new Font("Verdana", 8.25F),
                BackColor = Color.White,
                ForeColor = Color.Black,
                BorderStyle = BorderStyle.FixedSingle,
                TabIndex = 10
            };

            visualizationPanel = new Panel
            {
                Location = new Point(318, 12),
                Size = new Size(560, 520),
                TabIndex = 11,
                BackColor = Color.White,
                BorderStyle = BorderStyle.FixedSingle
            };
            visualizationPanel.Paint += visualizationPanel_Paint;

            controlPanel.Controls.AddRange(new Control[]
            {
                riverDischargeLabel, riverDischargeTextBox, tidalAmplitudeLabel, tidalAmplitudeTextBox,
                tidalPeriodLabel, tidalPeriodTextBox, windSpeedLabel, windSpeedTextBox,
                windDirectionLabel, windDirectionTextBox, salinityGradientLabel, salinityGradientTextBox,
                waveHeightLabel, waveHeightTextBox, stormSurgeAmplitudeLabel, stormSurgeAmplitudeTextBox,
                solverTypeLabel, solverTypeComboBox, enableWetDryCheckBox, minDepthLabel, minDepthTextBox,
                startButton, pauseButton, resetButton
            });

            this.Controls.AddRange(new Control[] { controlPanel, outputConsoleTextBox, visualizationPanel });
        }

        private Label CreateLabel(string text, int x, int y)
        {
            return new Label
            {
                AutoSize = true,
                Location = new Point(x, y),
                Text = text,
                Font = new Font("Verdana", 8.25F),
                BackColor = SystemColors.Control,
                ForeColor = Color.Black
            };
        }

        private TextBox CreateTextBox(string text, int x, int y)
        {
            return new TextBox
            {
                Location = new Point(x, y),
                Size = new Size(250, 22),
                Text = text,
                Font = new Font("Verdana", 8.25F),
                BackColor = Color.White,
                ForeColor = Color.Black,
                BorderStyle = BorderStyle.FixedSingle
            };
        }

        private Button CreateButton(string text, int x, int y, EventHandler clickHandler)
        {
            var button = new Button
            {
                Location = new Point(x, y),
                Size = new Size(170, 25),
                Text = text,
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                BackColor = SystemColors.Control,
                ForeColor = Color.Black,
                Font = new Font("Verdana", 8.25F)
            };
            button.Click += clickHandler;
            return button;
        }

        private void InitializeSimulation()
        {
            solverType = solverTypeComboBox.SelectedItem?.ToString() ?? "1D Simplified";
            windForcing = new WindForcing(windSpeed, windDirection, waveHeight);

            /**
            if (solverType == "3D Navier-Stokes")
            {
                dx = estuaryLength / (nx - 1);
                dy = estuaryWidth / (ny - 1);
                dz = estuaryDepth / (nz - 1);
                u = new double[nx, ny, nz];
                v = new double[nx, ny, nz];
                w = new double[nx, ny, nz];
                p = new double[nx, ny, nz];
                salinity = new double[nx, ny, nz];
                temperature = new double[nx, ny, nz];

                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        for (int k = 0; k < nz; k++)
                        {
                            u[i, j, k] = riverDischarge / (estuaryWidth * estuaryDepth);
                            v[i, j, k] = 0.0;
                            w[i, j, k] = 0.0;
                            p[i, j, k] = 0.0;
                            salinity[i, j, k] = salinityGradient * i * dx;
                            temperature[i, j, k] = 20.0;
                        }
            }
            **/
            if (solverType == "3D Navier-Stokes")
            {
                dx = estuaryLength / (nx - 1);
                dy = estuaryWidth / (ny - 1);
                dz = estuaryDepth / (nz - 1);
                u = new double[nx, ny, nz];
                v = new double[nx, ny, nz];
                w = new double[nx, ny, nz];
                p = new double[nx, ny, nz];
                salinity = new double[nx, ny, nz];
                temperature = new double[nx, ny, nz];

                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        for (int k = 0; k < nz; k++)
                        {
                            u[i, j, k] = riverDischarge / (estuaryWidth * estuaryDepth);
                            v[i, j, k] = 0.0;
                            w[i, j, k] = 0.0;
                            // Initialize pressure to reflect tidal and surge effects
                            double x = i * dx;
                            double eta = tidalAmplitude * Math.Sin(2 * Math.PI * x / estuaryLength) + stormSurgeAmplitude * Math.Exp(-x / estuaryLength);
                            p[i, j, k] = rho0 * g * eta * (1.0 - k * dz / estuaryDepth); // Hydrostatic pressure
                            salinity[i, j, k] = salinityGradient * x;
                            temperature[i, j, k] = 20.0;
                        }
            }
            else if (solverType == "2D Shallow Water")
            {
                swe2DSolver = new ShallowWaterEq2D(estuaryLength, estuaryWidth, estuaryDepth, nx, ny, 
                    riverDischarge, tidalAmplitude, tidalPeriod, windForcing, salinityGradient,
                    waveHeight, stormSurgeAmplitude, enableWetDry, minDepth);
                if (enableWetDry)
                {
                    double[,] bathymetry = new double[nx, ny];
                    for (int i = 0; i < nx; i++)
                        for (int j = 0; j < ny; j++)
                        {
                            bathymetry[i, j] = estuaryDepth * (1.0 - i * dx / estuaryLength);
                        }
                    wetDryAlgo = new WetAndDryAlgo(bathymetry, minDepth, nx, ny, estuaryLength, estuaryWidth);
                    wetDryAlgo.InitializeWetDry(swe2DSolver.Eta);
                }
            }

            /*
            else // 1D Simplified
            {
                velocityProfile = new double[gridPoints];
                salinityProfile = new double[gridPoints];
                waterLevel = new double[gridPoints];
                double dx1D = estuaryLength / (gridPoints - 1);
                for (int i = 0; i < gridPoints; i++)
                {
                    velocityProfile[i] = riverDischarge / estuaryDepth;
                    salinityProfile[i] = salinityGradient * i * dx1D;
                    waterLevel[i] = tidalAmplitude * Math.Sin(2 * Math.PI * (i * dx1D) / estuaryLength);
                }
            }
            */
            else // 1D Simplified
            {
                velocityProfile = new double[gridPoints];
                salinityProfile = new double[gridPoints];
                waterLevel = new double[gridPoints];
                double dx1D = estuaryLength / (gridPoints - 1);
                for (int i = 0; i < gridPoints; i++)
                {
                    velocityProfile[i] = riverDischarge / estuaryDepth;
                    salinityProfile[i] = salinityGradient * i * dx1D;
                    waterLevel[i] = tidalAmplitude * Math.Sin(2 * Math.PI * (i * dx1D) / estuaryLength);
                }
                double initialAvgWaterLevel = waterLevel.Average();
                double initialMaxWaterLevel = waterLevel.Max(Math.Abs);
                outputConsoleTextBox.AppendText($"Initial 1D Water Level: Avg={initialAvgWaterLevel:F4} m, Max={initialMaxWaterLevel:F4} m\r\n");
            }
            if (simulationTimer != null)
            {
                simulationTimer.Stop();
                simulationTimer.Dispose();
            }
            simulationTimer = new Timer();
            simulationTimer.Interval = 100;
            simulationTimer.Tick += (s, e) => UpdateSimulation();
        }

        private double ComputeCFLTimeStep()
        {
            if (solverType == "3D Navier-Stokes")
            {
                double maxU = u.Cast<double>().Max(Math.Abs);
                double maxV = v.Cast<double>().Max(Math.Abs);
                double maxW = w.Cast<double>().Max(Math.Abs);
                double cflDt = courantNumber * Math.Min(dx / (maxU + 1e-10), Math.Min(dy / (maxV + 1e-10), dz / (maxW + 1e-10)));
                cflDt = Math.Min(cflDt, Math.Min(dx * dx / (2 * (nu + eddyViscosity)), Math.Min(dy * dy / (2 * (nu + eddyViscosity)), dz * dz / (2 * (nu + eddyViscosity)))));
                return Math.Max(0.01, Math.Min(dt, cflDt));
            }
            else if (solverType == "2D Shallow Water")
            {
                return swe2DSolver.ComputeCFLTimeStep();
            }
            else
            {
                double maxU = velocityProfile.Max(Math.Abs);
                double dx1D = estuaryLength / (gridPoints - 1);
                double cflDt = courantNumber * dx1D / (maxU + 1e-10);
                cflDt = Math.Min(cflDt, dx1D * dx1D / (2 * (nu + eddyViscosity)));
                return Math.Max(0.01, Math.Min(dt, cflDt));
            }
        }

        private double ApplyFluxLimiter(double q, double qUpwind, double qDownwind)
        {
            double r = (q - qUpwind + 1e-10) / (qDownwind - q + 1e-10);
            double phi = Math.Max(0, Math.Min(1, r)); // Minmod limiter
            return phi * (qDownwind - q);
        }

        private void UpdateSimulation()
        {
            if (IsDisposed || outputConsoleTextBox.IsDisposed || visualizationPanel.IsDisposed)
                return;

            dt = ComputeCFLTimeStep();
            currentTime += dt;
            // Update windForcing parameters
            windForcing.UpdateParameters(windSpeed, windDirection, waveHeight);
            
            if (solverType == "3D Navier-Stokes")
            {
                double[,,] uNew = new double[nx, ny, nz];
                double[,,] vNew = new double[nx, ny, nz];
                double[,,] wNew = new double[nx, ny, nz];
                double[,,] salinityNew = new double[nx, ny, nz];
                double[,,] temperatureNew = new double[nx, ny, nz];

                double tidalVelocity = tidalAmplitude * (2 * Math.PI / tidalPeriod) * Math.Cos(2 * Math.PI * currentTime / tidalPeriod);
                double tauX = 0.001 * rho0 * windSpeed * windSpeed * Math.Cos(windDirection * Math.PI / 180.0);
                double tauY = 0.001 * rho0 * windSpeed * windSpeed * Math.Sin(windDirection * Math.PI / 180.0);
                //var (tauX, tauY) = windForcing.ComputeWindStress(); // Use WindForcing for wind stress
                (tauX, tauY) = windForcing.ComputeWindStress(); // Compute wind stress once
                double surgeLevel = stormSurgeAmplitude * Math.Exp(-currentTime / 86400.0);
                double tideSurgeInteraction = 0.05 * tidalAmplitude * stormSurgeAmplitude * Math.Sin(2 * Math.PI * currentTime / tidalPeriod);
                double waveNumber = 2 * Math.PI / (estuaryLength / 10);
                double waveVelocityX = 0.5 * waveHeight * waveNumber * Math.Sqrt(g * estuaryDepth);

                for (int i = 1; i < nx - 1; i++)
                    for (int j = 1; j < ny - 1; j++)
                        for (int k = 1; k < nz - 1; k++)
                        {
                            double rho = rho0 + salinity[i, j, k] * (rhoOcean - rho0) / 35.0 - 0.2 * (temperature[i, j, k] - 20.0);
                            double uAdvX = u[i, j, k] >= 0 ? u[i, j, k] * (u[i, j, k] - u[i - 1, j, k]) / dx : u[i, j, k] * (u[i + 1, j, k] - u[i, j, k]) / dx;
                            double uAdvY = v[i, j, k] >= 0 ? v[i, j, k] * (u[i, j, k] - u[i, j - 1, k]) / dy : v[i, j, k] * (u[i, j + 1, k] - u[i, j, k]) / dy;
                            double uAdvZ = w[i, j, k] >= 0 ? w[i, j, k] * (u[i, j, k] - u[i, j, k - 1]) / dz : w[i, j, k] * (u[i, j, k + 1] - u[i, j, k]) / dz;
                            double vAdvX = u[i, j, k] >= 0 ? u[i, j, k] * (v[i, j, k] - v[i - 1, j, k]) / dx : u[i, j, k] * (v[i + 1, j, k] - v[i, j, k]) / dx;
                            double vAdvY = v[i, j, k] >= 0 ? v[i, j, k] * (v[i, j, k] - v[i, j - 1, k]) / dy : v[i, j, k] * (v[i, j + 1, k] - v[i, j, k]) / dy;
                            double vAdvZ = w[i, j, k] >= 0 ? w[i, j, k] * (v[i, j, k] - v[i, j, k - 1]) / dz : w[i, j, k] * (v[i, j, k + 1] - v[i, j, k]) / dz;
                            double wAdvX = u[i, j, k] >= 0 ? u[i, j, k] * (w[i, j, k] - w[i - 1, j, k]) / dx : u[i, j, k] * (w[i + 1, j, k] - w[i, j, k]) / dx;
                            double wAdvY = v[i, j, k] >= 0 ? v[i, j, k] * (w[i, j, k] - w[i, j - 1, k]) / dy : v[i, j, k] * (w[i, j + 1, k] - w[i, j, k]) / dy;
                            double wAdvZ = w[i, j, k] >= 0 ? w[i, j, k] * (w[i, j, k] - w[i, j, k - 1]) / dz : w[i, j, k] * (w[i, j, k + 1] - w[i, j, k]) / dz;

                            uAdvX += u[i, j, k] >= 0 ? ApplyFluxLimiter(u[i, j, k], u[i - 1, j, k], u[i + 1, j, k]) / dx : ApplyFluxLimiter(u[i, j, k], u[i + 1, j, k], u[i - 1, j, k]) / dx;
                            uAdvY += v[i, j, k] >= 0 ? ApplyFluxLimiter(u[i, j, k], u[i, j - 1, k], u[i, j + 1, k]) / dy : ApplyFluxLimiter(u[i, j, k], u[i, j + 1, k], u[i, j - 1, k]) / dy;
                            uAdvZ += w[i, j, k] >= 0 ? ApplyFluxLimiter(u[i, j, k], u[i, j, k - 1], u[i, j, k + 1]) / dz : ApplyFluxLimiter(u[i, j, k], u[i, j, k + 1], u[i, j, k - 1]) / dz;
                            vAdvX += u[i, j, k] >= 0 ? ApplyFluxLimiter(v[i, j, k], v[i - 1, j, k], v[i + 1, j, k]) / dx : ApplyFluxLimiter(v[i, j, k], v[i + 1, j, k], v[i - 1, j, k]) / dx;
                            vAdvY += v[i, j, k] >= 0 ? ApplyFluxLimiter(v[i, j, k], v[i, j - 1, k], v[i, j + 1, k]) / dy : ApplyFluxLimiter(v[i, j, k], v[i, j + 1, k], v[i, j - 1, k]) / dy;
                            vAdvZ += w[i, j, k] >= 0 ? ApplyFluxLimiter(v[i, j, k], v[i, j, k - 1], v[i, j, k + 1]) / dz : ApplyFluxLimiter(v[i, j, k], v[i, j, k + 1], v[i, j, k - 1]) / dz;
                            wAdvX += u[i, j, k] >= 0 ? ApplyFluxLimiter(w[i, j, k], w[i - 1, j, k], w[i + 1, j, k]) / dx : ApplyFluxLimiter(w[i, j, k], w[i + 1, j, k], w[i - 1, j, k]) / dx;
                            wAdvY += v[i, j, k] >= 0 ? ApplyFluxLimiter(w[i, j, k], w[i, j - 1, k], w[i, j + 1, k]) / dy : ApplyFluxLimiter(w[i, j, k], w[i, j + 1, k], w[i, j - 1, k]) / dy;
                            wAdvZ += w[i, j, k] >= 0 ? ApplyFluxLimiter(w[i, j, k], w[i, j, k - 1], w[i, j, k + 1]) / dz : ApplyFluxLimiter(w[i, j, k], w[i, j, k + 1], w[i, j, k - 1]) / dz;

                            double uDiff = (nu + eddyViscosity) * ((u[i + 1, j, k] - 2 * u[i, j, k] + u[i - 1, j, k]) / (dx * dx) +
                                                                   (u[i, j + 1, k] - 2 * u[i, j, k] + u[i, j - 1, k]) / (dy * dy) +
                                                                   (u[i, j, k + 1] - 2 * u[i, j, k] + u[i, j, k - 1]) / (dz * dz));
                            double vDiff = (nu + eddyViscosity) * ((v[i + 1, j, k] - 2 * v[i, j, k] + v[i - 1, j, k]) / (dx * dx) +
                                                                   (v[i, j + 1, k] - 2 * v[i, j, k] + v[i, j - 1, k]) / (dy * dy) +
                                                                   (v[i, j, k + 1] - 2 * v[i, j, k] + v[i, j, k - 1]) / (dz * dz));
                            double wDiff = (nu + eddyViscosity) * ((w[i + 1, j, k] - 2 * w[i, j, k] + w[i - 1, j, k]) / (dx * dx) +
                                                                   (w[i, j + 1, k] - 2 * w[i, j, k] + w[i, j - 1, k]) / (dy * dy) +
                                                                   (w[i, j, k + 1] - 2 * w[i, j, k] + w[i, j, k - 1]) / (dz * dz));
                            double dpdx = (p[i + 1, j, k] - p[i - 1, j, k]) / (2 * dx) + atmosphericPressureGradient / rho0;
                            double dpdy = (p[i, j + 1, k] - p[i, j - 1, k]) / (2 * dy);
                            double dpdz = (p[i, j, k + 1] - p[i, j, k - 1]) / (2 * dz);
                            double drhodx = i < nx - 1 && i > 0 ? (rho0 + salinity[i + 1, j, k] * (rhoOcean - rho0) / 35.0 - 0.2 * (temperature[i + 1, j, k] - 20.0) -
                                                                   (rho0 + salinity[i - 1, j, k] * (rhoOcean - rho0) / 35.0 - 0.2 * (temperature[i - 1, j, k] - 20.0))) / (2 * dx) : 0;
                            double baroclinicX = -g * drhodx / rho0;
                            double coriolisU = coriolisParameter * v[i, j, k];
                            double coriolisV = -coriolisParameter * u[i, j, k];
                            double quadraticDrag = -frictionCoefficient * Math.Abs(u[i, j, k]) * u[i, j, k] / estuaryDepth;
                            uNew[i, j, k] = u[i, j, k] + dt * (-uAdvX - uAdvY - uAdvZ - dpdx / rho0 + uDiff + baroclinicX + coriolisU + waveVelocityX + quadraticDrag);
                            vNew[i, j, k] = v[i, j, k] + dt * (-vAdvX - vAdvY - vAdvZ - dpdy / rho0 + vDiff + coriolisV);
                            wNew[i, j, k] = w[i, j, k] + dt * (-wAdvX - wAdvY - wAdvZ - dpdz / rho0 + wDiff - g);
                            if (k == nz - 1)
                            {
                                uNew[i, j, k] += dt * tauX / (rho0 * dz);
                                vNew[i, j, k] += dt * tauY / (rho0 * dz);
                            }
                        }

                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++)
                    {
                        uNew[0, j, k] = riverDischarge / (estuaryWidth * estuaryDepth);
                        vNew[0, j, k] = 0.0;
                        wNew[0, j, k] = 0.0;
                        salinityNew[0, j, k] = 0.0;
                        temperatureNew[0, j, k] = 20.0;
                        uNew[nx - 1, j, k] = tidalVelocity;
                        vNew[nx - 1, j, k] = 0.0;
                        wNew[nx - 1, j, k] = 0.0;
                        p[nx - 1, j, k] = rho0 * g * (surgeLevel + tideSurgeInteraction);
                        salinityNew[nx - 1, j, k] = 35.0 + seasonalSalinityAmplitude * Math.Sin(2 * Math.PI * currentTime / (365 * 86400));
                        temperatureNew[nx - 1, j, k] = 20.0;
                    }

                for (int i = 0; i < nx; i++)
                    for (int k = 0; k < nz; k++)
                    {
                        vNew[i, 0, k] = 0.0;
                        vNew[i, ny - 1, k] = 0.0;
                        uNew[i, 0, k] = u[i, 1, k];
                        uNew[i, ny - 1, k] = u[i, ny - 2, k];
                        wNew[i, 0, k] = w[i, 1, k];
                        wNew[i, ny - 1, k] = w[i, ny - 2, k];
                        salinityNew[i, 0, k] = salinity[i, 1, k];
                        salinityNew[i, ny - 1, k] = salinity[i, ny - 2, k];
                        temperatureNew[i, 0, k] = temperature[i, 1, k];
                        temperatureNew[i, ny - 1, k] = temperature[i, ny - 2, k];
                    }

                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        uNew[i, j, 0] = 0.0;
                        vNew[i, j, 0] = 0.0;
                        wNew[i, j, 0] = 0.0;
                    }

                double[,,] divergence = new double[nx, ny, nz];
                double maxDivergence = 0.0;
                for (int i = 1; i < nx - 1; i++)
                    for (int j = 1; j < ny - 1; j++)
                        for (int k = 1; k < nz - 1; k++)
                        {
                            divergence[i, j, k] = (uNew[i + 1, j, k] - uNew[i - 1, j, k]) / (2 * dx) +
                                                  (vNew[i, j + 1, k] - vNew[i, j - 1, k]) / (2 * dy) +
                                                  (wNew[i, j, k + 1] - wNew[i, j, k - 1]) / (2 * dz);
                            maxDivergence = Math.Max(maxDivergence, Math.Abs(divergence[i, j, k]));
                        }

                for (int iter = 0; iter < 500; iter++)
                {
                    double maxChange = 0.0;
                    for (int i = 1; i < nx - 1; i++)
                        for (int j = 1; j < ny - 1; j++)
                            for (int k = 1; k < nz - 1; k++)
                            {
                                double pOld = p[i, j, k];
                                p[i, j, k] = (p[i + 1, j, k] + p[i - 1, j, k] + p[i, j + 1, k] + p[i, j - 1, k] + p[i, j, k + 1] + p[i, j, k - 1] -
                                              dx * dx * divergence[i, j, k]) / 6.0;
                                maxChange = Math.Max(maxChange, Math.Abs(p[i, j, k] - pOld));
                            }
                    if (maxChange < 1e-6) break;
                }

                for (int i = 1; i < nx - 1; i++)
                    for (int j = 1; j < ny - 1; j++)
                        for (int k = 1; k < nz - 1; k++)
                        {
                            uNew[i, j, k] -= dt * (p[i + 1, j, k] - p[i - 1, j, k]) / (2 * dx * rho0);
                            vNew[i, j, k] -= dt * (p[i, j + 1, k] - p[i, j - 1, k]) / (2 * dy * rho0);
                            wNew[i, j, k] -= dt * (p[i, j, k + 1] - p[i, j, k - 1]) / (2 * dz * rho0);

                            double sAdvX = u[i, j, k] >= 0 ? u[i, j, k] * (salinity[i, j, k] - salinity[i - 1, j, k]) / dx : u[i, j, k] * (salinity[i + 1, j, k] - salinity[i, j, k]) / dx;
                            double sAdvY = v[i, j, k] >= 0 ? v[i, j, k] * (salinity[i, j, k] - salinity[i, j - 1, k]) / dy : v[i, j, k] * (salinity[i, j + 1, k] - salinity[i, j, k]) / dy;
                            double sAdvZ = w[i, j, k] >= 0 ? w[i, j, k] * (salinity[i, j, k] - salinity[i, j, k - 1]) / dz : w[i, j, k] * (salinity[i, j, k + 1] - salinity[i, j, k]) / dz;
                            sAdvX += u[i, j, k] >= 0 ? ApplyFluxLimiter(salinity[i, j, k], salinity[i - 1, j, k], salinity[i + 1, j, k]) / dx : ApplyFluxLimiter(salinity[i, j, k], salinity[i + 1, j, k], salinity[i - 1, j, k]) / dx;
                            sAdvY += v[i, j, k] >= 0 ? ApplyFluxLimiter(salinity[i, j, k], salinity[i, j - 1, k], salinity[i, j + 1, k]) / dy : ApplyFluxLimiter(salinity[i, j, k], salinity[i, j + 1, k], salinity[i, j - 1, k]) / dy;
                            sAdvZ += w[i, j, k] >= 0 ? ApplyFluxLimiter(salinity[i, j, k], salinity[i, j, k - 1], salinity[i, j, k + 1]) / dz : ApplyFluxLimiter(salinity[i, j, k], salinity[i, j, k + 1], salinity[i, j, k - 1]) / dz;
                            double sDiff = kappa * ((salinity[i + 1, j, k] - 2 * salinity[i, j, k] + salinity[i - 1, j, k]) / (dx * dx) +
                                                    (salinity[i, j + 1, k] - 2 * salinity[i, j, k] + salinity[i, j - 1, k]) / (dy * dy) +
                                                    (salinity[i, j, k + 1] - 2 * salinity[i, j, k] + salinity[i, j, k - 1]) / (dz * dz));
                            double shearMixing = 0.1 * Math.Abs((u[i + 1, j, k] - u[i - 1, j, k]) / (2 * dx)) * salinity[i, j, k];
                            salinityNew[i, j, k] = salinity[i, j, k] + dt * (-sAdvX - sAdvY - sAdvZ + sDiff + shearMixing);
                            salinityNew[i, j, k] = Math.Max(0.0, Math.Min(35.0, salinityNew[i, j, k]));

                            double tAdvX = u[i, j, k] >= 0 ? u[i, j, k] * (temperature[i, j, k] - temperature[i - 1, j, k]) / dx : u[i, j, k] * (temperature[i + 1, j, k] - temperature[i, j, k]) / dx;
                            double tAdvY = v[i, j, k] >= 0 ? v[i, j, k] * (temperature[i, j, k] - temperature[i, j - 1, k]) / dy : v[i, j, k] * (temperature[i, j + 1, k] - temperature[i, j, k]) / dy;
                            double tAdvZ = w[i, j, k] >= 0 ? w[i, j, k] * (temperature[i, j, k] - temperature[i, j, k - 1]) / dz : w[i, j, k] * (temperature[i, j, k + 1] - temperature[i, j, k]) / dz;
                            tAdvX += u[i, j, k] >= 0 ? ApplyFluxLimiter(temperature[i, j, k], temperature[i - 1, j, k], temperature[i + 1, j, k]) / dx : ApplyFluxLimiter(temperature[i, j, k], temperature[i + 1, j, k], temperature[i - 1, j, k]) / dx;
                            tAdvY += v[i, j, k] >= 0 ? ApplyFluxLimiter(temperature[i, j, k], temperature[i, j - 1, k], temperature[i, j + 1, k]) / dy : ApplyFluxLimiter(temperature[i, j, k], temperature[i, j + 1, k], temperature[i, j - 1, k]) / dy;
                            tAdvZ += w[i, j, k] >= 0 ? ApplyFluxLimiter(temperature[i, j, k], temperature[i, j, k - 1], temperature[i, j, k + 1]) / dz : ApplyFluxLimiter(temperature[i, j, k], temperature[i, j, k + 1], temperature[i, j, k - 1]) / dz;
                            double tDiff = kappa * ((temperature[i + 1, j, k] - 2 * temperature[i, j, k] + temperature[i - 1, j, k]) / (dx * dx) +
                                                    (temperature[i, j + 1, k] - 2 * temperature[i, j, k] + temperature[i, j - 1, k]) / (dy * dy) +
                                                    (temperature[i, j, k + 1] - 2 * temperature[i, j, k] + temperature[i, j, k - 1]) / (dz * dz));
                            temperatureNew[i, j, k] = temperature[i, j, k] + dt * (-tAdvX - tAdvY - tAdvZ + tDiff);
                            temperatureNew[i, j, k] = Math.Max(0.0, Math.Min(40.0, temperatureNew[i, j, k]));
                        }

                u = uNew;
                v = vNew;
                w = wNew;
                salinity = salinityNew;
                temperature = temperatureNew;

                double maxU = u.Cast<double>().Max(Math.Abs);
                double maxV = v.Cast<double>().Max(Math.Abs);
                double maxW = w.Cast<double>().Max(Math.Abs);
                if (!outputConsoleTextBox.IsDisposed)
                    outputConsoleTextBox.AppendText($"CFL dt: {dt:F4}s | Max |u|: {maxU:F4} m/s | Max |v|: {maxV:F4} m/s | Max |w|: {maxW:F4} m/s | Max Div: {maxDivergence:F4}\r\n");

                if (u.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)) ||
                    v.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)) ||
                    w.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)) ||
                    salinity.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)) ||
                    temperature.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)))
                {
                    if (!outputConsoleTextBox.IsDisposed)
                        outputConsoleTextBox.AppendText("Warning: Numerical instability detected in 3D solver!\r\n");
                    simulationTimer.Stop();
                    isSimulationRunning = false;
                    UpdateButtonStates();
                    return;
                }
            }
            else if (solverType == "2D Shallow Water")
            {
                swe2DSolver.Update(currentTime, wetDryAlgo);
                if (enableWetDry && wetDryAlgo != null)
                {
                    wetDryAlgo.ApplyWetDry(swe2DSolver.U, swe2DSolver.V, swe2DSolver.Eta, swe2DSolver.Salinity, dt);
                }
                dt = swe2DSolver.ComputeCFLTimeStep();
                double maxU = swe2DSolver.U.Cast<double>().Max(Math.Abs);
                double maxV = swe2DSolver.V.Cast<double>().Max(Math.Abs);
                double minSalinity = swe2DSolver.Salinity.Cast<double>().Min();
                double maxSalinity = swe2DSolver.Salinity.Cast<double>().Max();
                double maxEta = swe2DSolver.Eta.Cast<double>().Max(Math.Abs);
                if (!outputConsoleTextBox.IsDisposed)
                    outputConsoleTextBox.AppendText($"CFL dt: {dt:F4}s | Max |u|: {maxU:F4} m/s | Max |v|: {maxV:F4} m/s | Min Salinity: {minSalinity:F2} PSU | Max Salinity: {maxSalinity:F2} PSU | Max Water Level: {maxEta:F2} m\r\n");

                if (swe2DSolver.U.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)) ||
                    swe2DSolver.V.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)) ||
                    swe2DSolver.Eta.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)) ||
                    swe2DSolver.Salinity.Cast<double>().Any(val => double.IsNaN(val) || double.IsInfinity(val)))
                {
                    if (!outputConsoleTextBox.IsDisposed)
                        outputConsoleTextBox.AppendText("Warning: Numerical instability detected in 2D solver!\r\n");
                    simulationTimer.Stop();
                    isSimulationRunning = false;
                    UpdateButtonStates();
                    return;
                }
            }
            // 1D Simplified Solver

            else
            {
                double dx1D = estuaryLength / (gridPoints - 1);
                double[] newVelocity = new double[gridPoints];
                double[] newWaterLevel = new double[gridPoints];
                double[] newSalinity = new double[gridPoints];

                double tidalForcing = tidalAmplitude * Math.Sin(2 * Math.PI * currentTime / tidalPeriod);
                double tidalVelocity = tidalAmplitude * (2 * Math.PI / tidalPeriod) * Math.Cos(2 * Math.PI * currentTime / tidalPeriod);
                var (tauX, _) = windForcing.ComputeWindStress(); // Use WindForcing for wind stress (x-component only for 1D)
                double windStress = 0.001 * rho0 * windSpeed * windSpeed * Math.Cos(windDirection * Math.PI / 180.0);
                double waveNumber = 2 * Math.PI / (estuaryLength / 10);
                double waveVelocity = 0.5 * waveHeight * waveNumber * Math.Sqrt(g * estuaryDepth);
                double surgeLevel = stormSurgeAmplitude * Math.Exp(-currentTime / 86400.0);
                double tideSurgeInteraction = 0.05 * tidalAmplitude * stormSurgeAmplitude * Math.Sin(2 * Math.PI * currentTime / tidalPeriod);

                newSalinity[0] = 0.0;
                newSalinity[gridPoints - 1] = 35.0 + seasonalSalinityAmplitude * Math.Sin(2 * Math.PI * currentTime / (365 * 86400));
                newVelocity[0] = riverDischarge / estuaryDepth;
                newVelocity[gridPoints - 1] = tidalVelocity;
                newWaterLevel[gridPoints - 1] = tidalForcing + surgeLevel;

                for (int i = 1; i < gridPoints - 1; i++)
                {
                    double rho = rho0 + salinityProfile[i] * (rhoOcean - rho0) / 35.0;
                    double drhodx = i < gridPoints - 1 && i > 0 ? (salinityProfile[i + 1] - salinityProfile[i - 1]) / (2 * dx1D) * (rhoOcean - rho0) / 35.0 : 0;
                    double dhdx = i < gridPoints - 1 && i > 0 ? (waterLevel[i + 1] - waterLevel[i - 1]) / (2 * dx1D) : 0;
                    double dudx = i < gridPoints - 1 && i > 0 ? (velocityProfile[i + 1] - velocityProfile[i - 1]) / (2 * dx1D) : 0;
                    double baroclinic = -g * drhodx / rho0;
                    double quadraticDrag = -frictionCoefficient * Math.Abs(velocityProfile[i]) * velocityProfile[i] / estuaryDepth;
                    double coriolis = coriolisParameter * velocityProfile[i]; // Simplified for 1D
                    double pressureGradient = atmosphericPressureGradient / rho0;
                    double shearMixing = 0.1 * Math.Abs(dudx) * salinityProfile[i];
                    double forcing = -g * dhdx + baroclinic + windStress / (rho0 * estuaryDepth) + waveVelocity + tideSurgeInteraction + quadraticDrag + coriolis + pressureGradient;
                    double uAdv = velocityProfile[i] >= 0 ? velocityProfile[i] * (velocityProfile[i] - velocityProfile[i - 1]) / dx1D : velocityProfile[i] * (velocityProfile[i + 1] - velocityProfile[i]) / dx1D;
                    uAdv += velocityProfile[i] >= 0 ? ApplyFluxLimiter(velocityProfile[i], velocityProfile[i - 1], velocityProfile[i + 1]) / dx1D : ApplyFluxLimiter(velocityProfile[i], velocityProfile[i + 1], velocityProfile[i - 1]) / dx1D;
                    newVelocity[i] = velocityProfile[i] + dt * (-uAdv + forcing);
                    newVelocity[i] = Math.Max(-1.0, Math.Min(1.0, newVelocity[i]));

                    double tidalWave = tidalAmplitude * Math.Sin(2 * Math.PI * (currentTime / tidalPeriod - i * dx1D / estuaryLength));
                    newWaterLevel[i] = waterLevel[i] - dt * (estuaryDepth + waterLevel[i]) * dudx + dt * tidalWave - dt * waterLevelDamping * waterLevel[i];

                    double dsdx = i < gridPoints - 1 && i > 0 ? (salinityProfile[i + 1] - salinityProfile[i - 1]) / (2 * dx1D) : 0;
                    double sAdv = velocityProfile[i] >= 0 ? velocityProfile[i] * (salinityProfile[i] - salinityProfile[i - 1]) / dx1D : velocityProfile[i] * (salinityProfile[i + 1] - salinityProfile[i]) / dx1D;
                    sAdv += velocityProfile[i] >= 0 ? ApplyFluxLimiter(salinityProfile[i], salinityProfile[i - 1], salinityProfile[i + 1]) / dx1D : ApplyFluxLimiter(salinityProfile[i], salinityProfile[i + 1], salinityProfile[i - 1]) / dx1D;
                    double sDiff = (kappa + eddyViscosity) * (salinityProfile[i + 1] - 2 * salinityProfile[i] + salinityProfile[i - 1]) / (dx1D * dx1D);
                    newSalinity[i] = salinityProfile[i] + dt * (-sAdv + sDiff + shearMixing);
                    newSalinity[i] = Math.Max(0.0, Math.Min(35.0, newSalinity[i]));
                }

                newSalinity[1] = 0.75 * newSalinity[1] + 0.25 * newSalinity[0];
                newSalinity[gridPoints - 2] = 0.75 * newSalinity[gridPoints - 2] + 0.25 * newSalinity[gridPoints - 1];

                velocityProfile = newVelocity;
                waterLevel = newWaterLevel;
                salinityProfile = newSalinity;

                double maxU = velocityProfile.Max(Math.Abs);
                double minSalinity = salinityProfile.Min();
                double maxSalinity = salinityProfile.Max();
                double maxWaterLevel = waterLevel.Max(Math.Abs);
                double avgWaterLevel = waterLevel.Average();
                if (!outputConsoleTextBox.IsDisposed)
                    outputConsoleTextBox.AppendText($"CFL dt: {dt:F4}s | Max |u|: {maxU:F4} m/s | Min Salinity: {minSalinity:F2} PSU | Max Salinity: {maxSalinity:F2} PSU | Max Water Level: {maxWaterLevel:F2} m\r\n");

                if (velocityProfile.Any(s => double.IsNaN(s) || double.IsInfinity(s)) ||
                    salinityProfile.Any(s => double.IsNaN(s) || double.IsInfinity(s)) ||
                    waterLevel.Any(s => double.IsNaN(s) || double.IsInfinity(s)))
                {
                    if (!outputConsoleTextBox.IsDisposed)
                        outputConsoleTextBox.AppendText("Warning: Numerical instability detected in 1D solver!\r\n");
                    simulationTimer.Stop();
                    isSimulationRunning = false;
                    UpdateButtonStates();
                    return;
                }
            }


            /**
            else // 1D Simplified Solver
            {
                double dx1D = estuaryLength / (gridPoints - 1);
                double[] newVelocity = new double[gridPoints];
                double[] newWaterLevel = new double[gridPoints];
                double[] newSalinity = new double[gridPoints];

                double tidalForcing = tidalAmplitude * Math.Sin(2 * Math.PI * currentTime / tidalPeriod);
                double tidalVelocity = tidalAmplitude * (2 * Math.PI / tidalPeriod) * Math.Cos(2 * Math.PI * currentTime / tidalPeriod);
                double windStress = 0.001 * rho0 * windSpeed * windSpeed * Math.Cos(windDirection * Math.PI / 180.0);
                double waveNumber = 2 * Math.PI / (estuaryLength / 10);
                double waveVelocity = 0.5 * waveHeight * waveNumber * Math.Sqrt(g * estuaryDepth);
                double surgeLevel = stormSurgeAmplitude * Math.Exp(-currentTime / 86400.0);
                double tideSurgeInteraction = 0.05 * tidalAmplitude * stormSurgeAmplitude * Math.Sin(2 * Math.PI * currentTime / tidalPeriod);

                newSalinity[0] = 0.0;
                newSalinity[gridPoints - 1] = 35.0 + seasonalSalinityAmplitude * Math.Sin(2 * Math.PI * currentTime / (365 * 86400));
                newVelocity[0] = riverDischarge / estuaryDepth;
                newVelocity[gridPoints - 1] = tidalVelocity;
                newWaterLevel[gridPoints - 1] = tidalForcing + surgeLevel + tideSurgeInteraction;

                for (int i = 1; i < gridPoints - 1; i++)
                {
                    double rho = rho0 + salinityProfile[i] * (rhoOcean - rho0) / 35.0;
                    double drhodx = i < gridPoints - 1 && i > 0 ? (salinityProfile[i + 1] - salinityProfile[i - 1]) / (2 * dx1D) * (rhoOcean - rho0) / 35.0 : 0;
                    double dhdx = i < gridPoints - 1 && i > 0 ? (waterLevel[i + 1] - waterLevel[i - 1]) / (2 * dx1D) : 0;
                    double dudx = i < gridPoints - 1 && i > 0 ? (velocityProfile[i + 1] - velocityProfile[i - 1]) / (2 * dx1D) : 0;
                    double baroclinic = -g * drhodx / rho0;
                    double quadraticDrag = -frictionCoefficient * Math.Abs(velocityProfile[i]) * velocityProfile[i] / estuaryDepth;
                    double coriolis = coriolisParameter * velocityProfile[i]; // Simplified for 1D
                    double pressureGradient = atmosphericPressureGradient / rho0;
                    double shearMixing = 0.1 * Math.Abs(dudx) * salinityProfile[i];
                    double forcing = -g * dhdx + baroclinic + windStress / (rho0 * estuaryDepth) + waveVelocity + tideSurgeInteraction + quadraticDrag + coriolis + pressureGradient;
                    double uAdv = velocityProfile[i] >= 0 ? velocityProfile[i] * (velocityProfile[i] - velocityProfile[i - 1]) / dx1D : velocityProfile[i] * (velocityProfile[i + 1] - velocityProfile[i]) / dx1D;
                    uAdv += velocityProfile[i] >= 0 ? ApplyFluxLimiter(velocityProfile[i], velocityProfile[i - 1], velocityProfile[i + 1]) / dx1D : ApplyFluxLimiter(velocityProfile[i], velocityProfile[i + 1], velocityProfile[i - 1]) / dx1D;
                    newVelocity[i] = velocityProfile[i] + dt * (-uAdv + forcing);
                    newVelocity[i] = Math.Max(-1.0, Math.Min(1.0, newVelocity[i]));

                    double tidalWave = tidalAmplitude * Math.Sin(2 * Math.PI * (currentTime / tidalPeriod - i * dx1D / (Math.Sqrt(g * estuaryDepth))));
                    double continuity = -(estuaryDepth + waterLevel[i]) * dudx;
                    newWaterLevel[i] = waterLevel[i] + dt * (continuity + tidalWave - waterLevelDamping * waterLevel[i]);
                    newWaterLevel[i] = Math.Max(-estuaryDepth, Math.Min(estuaryDepth, newWaterLevel[i])); // Limit water level

                    double dsdx = i < gridPoints - 1 && i > 0 ? (salinityProfile[i + 1] - salinityProfile[i - 1]) / (2 * dx1D) : 0;
                    double sAdv = velocityProfile[i] >= 0 ? velocityProfile[i] * (salinityProfile[i] - salinityProfile[i - 1]) / dx1D : velocityProfile[i] * (salinityProfile[i + 1] - salinityProfile[i]) / dx1D;
                    sAdv += velocityProfile[i] >= 0 ? ApplyFluxLimiter(salinityProfile[i], salinityProfile[i - 1], salinityProfile[i + 1]) / dx1D : ApplyFluxLimiter(salinityProfile[i], salinityProfile[i + 1], salinityProfile[i - 1]) / dx1D;
                    double sDiff = (kappa + eddyViscosity) * (salinityProfile[i + 1] - 2 * salinityProfile[i] + salinityProfile[i - 1]) / (dx1D * dx1D);
                    newSalinity[i] = salinityProfile[i] + dt * (-sAdv + sDiff + shearMixing);
                    newSalinity[i] = Math.Max(0.0, Math.Min(35.0, newSalinity[i]));
                }

                newSalinity[1] = 0.75 * newSalinity[1] + 0.25 * newSalinity[0];
                newSalinity[gridPoints - 2] = 0.75 * newSalinity[gridPoints - 2] + 0.25 * newSalinity[gridPoints - 1];

                velocityProfile = newVelocity;
                waterLevel = newWaterLevel;
                salinityProfile = newSalinity;

                double maxU = velocityProfile.Max(Math.Abs);
                double minSalinity = salinityProfile.Min();
                double maxSalinity = salinityProfile.Max();
                double maxWaterLevel = waterLevel.Max(Math.Abs);
                double avgWaterLevel = waterLevel.Average();
                if (!outputConsoleTextBox.IsDisposed)
                    outputConsoleTextBox.AppendText($"CFL dt: {dt:F4}s | Max |u|: {maxU:F4} m/s | Min Salinity: {minSalinity:F2} PSU | Max Salinity: {maxSalinity:F2} PSU | Max Water Level: {maxWaterLevel:F2} m | Avg Water Level: {avgWaterLevel:F2} m\r\n");

                if (velocityProfile.Any(s => double.IsNaN(s) || double.IsInfinity(s)) ||
                    salinityProfile.Any(s => double.IsNaN(s) || double.IsInfinity(s)) ||
                    waterLevel.Any(s => double.IsNaN(s) || double.IsInfinity(s)))
                {
                    if (!outputConsoleTextBox.IsDisposed)
                        outputConsoleTextBox.AppendText("Warning: Numerical instability detected in 1D solver!\r\n");
                    simulationTimer.Stop();
                    isSimulationRunning = false;
                    UpdateButtonStates();
                    return;
                }
            }
            **/

            if (!visualizationPanel.IsDisposed)
                visualizationPanel.Invalidate();
            if (!outputConsoleTextBox.IsDisposed)
                UpdateOutputConsole();
        }

        private void startButton_Click(object sender, EventArgs e)
        {
            try
            {
                riverDischarge = Math.Max(0.01, Math.Min(100.0, double.Parse(riverDischargeTextBox.Text)));
                tidalAmplitude = Math.Max(0.1, Math.Min(10.0, double.Parse(tidalAmplitudeTextBox.Text)));
                tidalPeriod = Math.Max(3600.0, Math.Min(86400.0, double.Parse(tidalPeriodTextBox.Text)));
                windSpeed = Math.Max(0.0, Math.Min(50.0, double.Parse(windSpeedTextBox.Text)));
                windDirection = Math.Max(0.0, Math.Min(360.0, double.Parse(windDirectionTextBox.Text)));
                salinityGradient = Math.Max(0.0, Math.Min(0.01, double.Parse(salinityGradientTextBox.Text)));
                waveHeight = Math.Max(0.0, Math.Min(5.0, double.Parse(waveHeightTextBox.Text)));
                stormSurgeAmplitude = Math.Max(0.0, Math.Min(5.0, double.Parse(stormSurgeAmplitudeTextBox.Text)));
                minDepth = Math.Max(0.001, Math.Min(1.0, double.Parse(minDepthTextBox.Text)));

                InitializeSimulation();
                simulationTimer.Start();
                isSimulationRunning = true;
                UpdateButtonStates();
                if (!outputConsoleTextBox.IsDisposed)
                    outputConsoleTextBox.AppendText($"Simulation started ({solverType} solver, Wet/Dry: {enableWetDry}).\r\n");
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Error: {ex.Message}", "Input Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private void pauseButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isSimulationRunning = false;
            UpdateButtonStates();
            if (!outputConsoleTextBox.IsDisposed)
                outputConsoleTextBox.AppendText("Simulation paused.\r\n");
        }

        private void resetButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isSimulationRunning = false;
            currentTime = 0.0;
            InitializeSimulation();
            if (!visualizationPanel.IsDisposed)
                visualizationPanel.Invalidate();
            if (!outputConsoleTextBox.IsDisposed)
            {
                outputConsoleTextBox.Clear();
                outputConsoleTextBox.AppendText("Simulation reset.\r\n");
            }
            UpdateButtonStates();
        }

        private void UpdateButtonStates()
        {
            if (!IsDisposed)
            {
                startButton.Enabled = !isSimulationRunning;
                pauseButton.Enabled = isSimulationRunning;
                resetButton.Enabled = true;
            }
        }

        private void UpdateOutputConsole()
        {
            /**
            if (solverType == "3D Navier-Stokes")
            {
                double avgU = 0.0, avgSalinity = 0.0, avgWaterLevel = 0.0, avgTemperature = 0.0;
                double minSalinity = double.MaxValue, maxSalinity = double.MinValue;
                double maxWaterLevel = 0.0;
                int count = 0;
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        for (int k = 0; k < nz; k++)
                        {
                            avgU += u[i, j, k];
                            avgSalinity += salinity[i, j, k];
                            //avgWaterLevel += Math.Abs(w[i, j, k]);
                            avgTemperature += temperature[i, j, k];
                            minSalinity = Math.Min(minSalinity, salinity[i, j, k]);
                            maxSalinity = Math.Max(maxSalinity, salinity[i, j, k]);
                            if (k == nz - 1) // Only use surface pressure for water level
                            {
                                double eta = p[i, j, k] / (rho0 * g); // Surface elevation from pressure
                                avgWaterLevel += Math.Abs(eta); // Use absolute value to avoid cancellation
                                maxWaterLevel = Math.Max(maxWaterLevel, Math.Abs(eta));
                                count++; // Count only surface points
                            }
                            maxWaterLevel = Math.Max(maxWaterLevel, Math.Abs(w[i, j, k]));
                            count++;
                        }
                avgU /= count;
                avgSalinity /= count;
                avgWaterLevel /= count;
                avgTemperature /= count;
                if (!outputConsoleTextBox.IsDisposed)
                {
                    string output = $"Time: {currentTime:F2}s | Avg Velocity X: {avgU:F4} m/s | Avg Salinity: {avgSalinity:F2} PSU | Min Salinity: {minSalinity:F2} PSU | Max Salinity: {maxSalinity:F2} PSU | Avg Water Level: {avgWaterLevel:F2} m | Max Water Level: {maxWaterLevel:F2} m | Avg Temperature: {avgTemperature:F2} °C";
                    outputConsoleTextBox.AppendText(output + "\r\n");
                }
            }
            **/
            if (solverType == "3D Navier-Stokes")
            {
                double avgU = 0.0, avgSalinity = 0.0, avgWaterLevel = 0.0, avgTemperature = 0.0;
                double minSalinity = double.MaxValue, maxSalinity = double.MinValue;
                double maxWaterLevel = 0.0;
                int count = 0;
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        // Compute surface elevation at k = nz-1
                        double eta = p[i, j, nz - 1] / (rho0 * g);
                        avgWaterLevel += Math.Abs(eta);
                        maxWaterLevel = Math.Max(maxWaterLevel, Math.Abs(eta));
                        count++;
                        for (int k = 0; k < nz; k++)
                        {
                            avgU += u[i, j, k];
                            avgSalinity += salinity[i, j, k];
                            avgTemperature += temperature[i, j, k];
                            minSalinity = Math.Min(minSalinity, salinity[i, j, k]);
                            maxSalinity = Math.Max(maxSalinity, salinity[i, j, k]);
                        }
                    }
                avgU /= (nx * ny * nz);
                avgSalinity /= (nx * ny * nz);
                avgTemperature /= (nx * ny * nz);
                avgWaterLevel /= count; // Average over surface points
                if (!outputConsoleTextBox.IsDisposed)
                {
                    string output = $"Time: {currentTime:F2}s | Avg Velocity X: {avgU:F4} m/s | Avg Salinity: {avgSalinity:F2} PSU | Min Salinity: {minSalinity:F2} PSU | Max Salinity: {maxSalinity:F2} PSU | Avg Water Level: {avgWaterLevel:F2} m | Max Water Level: {maxWaterLevel:F2} m | Avg Temperature: {avgTemperature:F2} °C";
                    outputConsoleTextBox.AppendText(output + "\r\n");
                }
            }
            else if (solverType == "2D Shallow Water")
            {
                double avgU = 0.0, avgV = 0.0, avgSalinity = 0.0, avgEta = 0.0;
                double minSalinity = double.MaxValue, maxSalinity = double.MinValue;
                double maxEta = 0.0;
                int wetCells = 0;
                int count = 0;
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        avgU += swe2DSolver.U[i, j];
                        avgV += swe2DSolver.V[i, j];
                        avgSalinity += swe2DSolver.Salinity[i, j];
                        avgEta += swe2DSolver.Eta[i, j];
                        minSalinity = Math.Min(minSalinity, swe2DSolver.Salinity[i, j]);
                        maxSalinity = Math.Max(maxSalinity, swe2DSolver.Salinity[i, j]);
                        maxEta = Math.Max(maxEta, Math.Abs(swe2DSolver.Eta[i, j]));
                        if (enableWetDry && wetDryAlgo != null && wetDryAlgo.GetWetDryStatus()[i, j])
                            wetCells++;
                        count++;
                    }
                avgU /= count;
                avgV /= count;
                avgSalinity /= count;
                avgEta /= count;
                if (!outputConsoleTextBox.IsDisposed)
                {
                    string output = $"Time: {currentTime:F2}s | Avg Velocity X: {avgU:F4} m/s | Avg Velocity Y: {avgV:F4} m/s | Avg Salinity: {avgSalinity:F2} PSU | Min Salinity: {minSalinity:F2} PSU | Max Salinity: {maxSalinity:F2} PSU | Avg Water Level: {avgEta:F2} m | Max Water Level: {maxEta:F2} m";
                    if (enableWetDry)
                        output += $" | Wet Cells: {wetCells}/{count}";
                    outputConsoleTextBox.AppendText(output + "\r\n");
                }
            }
            /**
            else
            {
                double avgVelocity = 0.0, avgSalinity = 0.0, avgWaterLevel = 0.0;
                double minSalinity = double.MaxValue, maxSalinity = double.MinValue;
                double maxWaterLevel = 0.0;
                for (int i = 0; i < gridPoints; i++)
                {
                    avgVelocity += velocityProfile[i];
                    avgSalinity += salinityProfile[i];
                    avgWaterLevel += waterLevel[i];
                    minSalinity = Math.Min(minSalinity, salinityProfile[i]);
                    maxSalinity = Math.Max(maxSalinity, salinityProfile[i]);
                    maxWaterLevel = Math.Max(maxWaterLevel, Math.Abs(waterLevel[i]));
                }
                avgVelocity /= gridPoints;
                avgSalinity /= gridPoints;
                avgWaterLevel /= gridPoints;
                if (!outputConsoleTextBox.IsDisposed)
                {
                    string output = $"Time: {currentTime:F2}s | Avg Velocity: {avgVelocity:F4} m/s | Avg Salinity: {avgSalinity:F2} PSU | Min Salinity: {minSalinity:F2} PSU | Max Salinity: {maxSalinity:F2} PSU | Avg Water Level: {avgWaterLevel:F2} m | Max Water Level: {maxWaterLevel:F2} m";
                    outputConsoleTextBox.AppendText(output + "\r\n");
                }
            }
            **/
            else // 1D Simplified
            {
                double avgVelocity = 0.0, avgSalinity = 0.0, avgWaterLevel = 0.0;
                double minSalinity = double.MaxValue, maxSalinity = double.MinValue;
                double maxWaterLevel = 0.0;
                for (int i = 0; i < gridPoints; i++)
                {
                    avgVelocity += velocityProfile[i];
                    avgSalinity += salinityProfile[i];
                    avgWaterLevel += Math.Abs(waterLevel[i]); // Use absolute value to avoid cancellation
                    minSalinity = Math.Min(minSalinity, salinityProfile[i]);
                    maxSalinity = Math.Max(maxSalinity, salinityProfile[i]);
                    maxWaterLevel = Math.Max(maxWaterLevel, Math.Abs(waterLevel[i]));
                }
                avgVelocity /= gridPoints;
                avgSalinity /= gridPoints;
                avgWaterLevel /= gridPoints;
                if (!outputConsoleTextBox.IsDisposed)
                {
                    string output = $"Time: {currentTime:F2}s | Avg Velocity: {avgVelocity:F4} m/s | Avg Salinity: {avgSalinity:F2} PSU | Min Salinity: {minSalinity:F2} PSU | Max Salinity: {maxSalinity:F2} PSU | Avg Water Level: {avgWaterLevel:F2} m | Max Water Level: {maxWaterLevel:F2} m";
                    outputConsoleTextBox.AppendText(output + "\r\n");
                }
            }
        }

        private void visualizationPanel_Paint(object sender, PaintEventArgs e)
        {
            if (visualizationPanel.IsDisposed)
                return;

            Graphics g = e.Graphics;
            g.Clear(Color.White);
            float panelWidth = visualizationPanel.Width;
            float panelHeight = visualizationPanel.Height;

            if (solverType == "3D Navier-Stokes")
            {
                double[] uAvg = new double[nx];
                double[] salinityAvg = new double[nx];
                double[] waterLevelAvg = new double[nx];
                for (int i = 0; i < nx; i++)
                {
                    double uSum = 0.0, sSum = 0.0, wSum = 0.0;
                    int kCount = 0;
                    for (int k = 0; k < nz; k++)
                        for (int j = 0; j < ny; j++)
                        {
                            uSum += u[i, j, k];
                            sSum += salinity[i, j, k];
                            wSum += w[i, j, k];
                            kCount++;
                        }
                    uAvg[i] = uSum / kCount;
                    salinityAvg[i] = sSum / kCount;
                    waterLevelAvg[i] = wSum / kCount;
                }

                using (Pen pen = new Pen(Color.Blue, 2))
                {
                    float maxU = (float)uAvg.Max(Math.Abs);
                    if (maxU == 0 || float.IsNaN(maxU) || float.IsInfinity(maxU)) maxU = 1.0f;
                    for (int i = 0; i < nx - 1; i++)
                    {
                        float x1 = i * panelWidth / (nx - 1);
                        float x2 = (i + 1) * panelWidth / (nx - 1);
                        float y1 = panelHeight / 2 - (float)(uAvg[i] / maxU) * panelHeight / 4;
                        float y2 = panelHeight / 2 - (float)(uAvg[i + 1] / maxU) * panelHeight / 4;
                        if (float.IsNaN(y1) || float.IsInfinity(y1) || float.IsNaN(y2) || float.IsInfinity(y2)) continue;
                        y1 = Math.Max(0, Math.Min(panelHeight, y1));
                        y2 = Math.Max(0, Math.Min(panelHeight, y2));
                        g.DrawLine(pen, x1, y1, x2, y2);
                    }
                }

                using (Pen pen = new Pen(Color.Red, 2))
                {
                    float maxSalinity = (float)salinityAvg.Max();
                    if (maxSalinity == 0 || float.IsNaN(maxSalinity) || float.IsInfinity(maxSalinity)) maxSalinity = 35.0f;
                    for (int i = 0; i < nx - 1; i++)
                    {
                        float x1 = i * panelWidth / (nx - 1);
                        float x2 = (i + 1) * panelWidth / (nx - 1);
                        float y1 = panelHeight - (float)(salinityAvg[i] / maxSalinity) * panelHeight / 4 - panelHeight / 2;
                        float y2 = panelHeight - (float)(salinityAvg[i + 1] / maxSalinity) * panelHeight / 4 - panelHeight / 2;
                        if (float.IsNaN(y1) || float.IsInfinity(y1) || float.IsNaN(y2) || float.IsInfinity(y2)) continue;
                        y1 = Math.Max(0, Math.Min(panelHeight, y1));
                        y2 = Math.Max(0, Math.Min(panelHeight, y2));
                        g.DrawLine(pen, x1, y1, x2, y2);
                    }
                }

                using (Pen pen = new Pen(Color.Green, 2))
                {
                    float maxWaterLevel = (float)waterLevelAvg.Max(Math.Abs);
                    if (maxWaterLevel == 0 || float.IsNaN(maxWaterLevel) || float.IsInfinity(maxWaterLevel)) maxWaterLevel = tidalAmplitude > 0 ? (float)tidalAmplitude : 1.0f;
                    for (int i = 0; i < nx - 1; i++)
                    {
                        float x1 = i * panelWidth / (nx - 1);
                        float x2 = (i + 1) * panelWidth / (nx - 1);
                        float y1 = panelHeight / 2 - (float)(waterLevelAvg[i] / maxWaterLevel) * panelHeight / 8;
                        float y2 = panelHeight / 2 - (float)(waterLevelAvg[i + 1] / maxWaterLevel) * panelHeight / 8;
                        if (float.IsNaN(y1) || float.IsInfinity(y1) || float.IsNaN(y2) || float.IsInfinity(y2)) continue;
                        y1 = Math.Max(0, Math.Min(panelHeight, y1));
                        y2 = Math.Max(0, Math.Min(panelHeight, y2));
                        g.DrawLine(pen, x1, y1, x2, y2);
                    }
                }
            }
            if (solverType == "2D Shallow Water")
            {
                double[] uAvg = new double[nx];
                double[] vAvg = new double[nx];
                double[] salinityAvg = new double[nx];
                double[] etaAvg = new double[nx];
                for (int i = 0; i < nx; i++)
                {
                    double uSum = 0.0, vSum = 0.0, sSum = 0.0, etaSum = 0.0;
                    int count = 0;
                    for (int j = 0; j < ny; j++)
                    {
                        uSum += swe2DSolver.U[i, j];
                        vSum += swe2DSolver.V[i, j];
                        sSum += swe2DSolver.Salinity[i, j];
                        etaSum += swe2DSolver.Eta[i, j];
                        count++;
                    }
                    uAvg[i] = uSum / count;
                    vAvg[i] = vSum / count;
                    salinityAvg[i] = sSum / count;
                    etaAvg[i] = etaSum / count;
                }

                // Draw wet/dry cells as background
                if (enableWetDry && wetDryAlgo != null)
                {
                    bool[,] wetDryStatus = wetDryAlgo.GetWetDryStatus();
                    float cellWidth = panelWidth / nx;
                    float cellHeight = panelHeight / ny;
                    for (int i = 0; i < nx; i++)
                        for (int j = 0; j < ny; j++)
                        {
                            float x = i * cellWidth;
                            float y = j * cellHeight;
                            Brush brush = wetDryStatus[i, j] ? new SolidBrush(Color.FromArgb(100, 135, 206, 235)) : new SolidBrush(Color.FromArgb(100, 210, 180, 140));
                            g.FillRectangle(brush, x, y, cellWidth, cellHeight);
                            brush.Dispose();
                        }
                }

                // Plot velocity (x-component)
                using (Pen pen = new Pen(Color.Blue, 2))
                {
                    float maxU = (float)uAvg.Max(Math.Abs);
                    if (maxU == 0 || float.IsNaN(maxU) || float.IsInfinity(maxU)) maxU = 1.0f;
                    for (int i = 0; i < nx - 1; i++)
                    {
                        float x1 = i * panelWidth / (nx - 1);
                        float x2 = (i + 1) * panelWidth / (nx - 1);
                        float y1 = panelHeight / 2 - (float)(uAvg[i] / maxU) * panelHeight / 4;
                        float y2 = panelHeight / 2 - (float)(uAvg[i + 1] / maxU) * panelHeight / 4;
                        if (float.IsNaN(y1) || float.IsInfinity(y1) || float.IsNaN(y2) || float.IsInfinity(y2)) continue;
                        y1 = Math.Max(0, Math.Min(panelHeight, y1));
                        y2 = Math.Max(0, Math.Min(panelHeight, y2));
                        g.DrawLine(pen, x1, y1, x2, y2);
                    }
                }

                // Plot salinity
                using (Pen pen = new Pen(Color.Red, 2))
                {
                    float maxSalinity = (float)salinityProfile.Max();
                    if (maxSalinity == 0 || float.IsNaN(maxSalinity) || float.IsInfinity(maxSalinity)) maxSalinity = 35.0f;
                    for (int i = 0; i < gridPoints - 1; i++)
                    {
                        float x1 = i * panelWidth / (gridPoints - 1);
                        float x2 = (i + 1) * panelWidth / (gridPoints - 1);
                        float y1 = panelHeight - (float)(salinityProfile[i] / maxSalinity) * panelHeight / 4 - panelHeight / 2;
                        float y2 = panelHeight - (float)(salinityProfile[i + 1] / maxSalinity) * panelHeight / 4 - panelHeight / 2;
                        if (float.IsNaN(y1) || float.IsInfinity(y1) || float.IsNaN(y2) || float.IsInfinity(y2)) continue;
                        y1 = Math.Max(0, Math.Min(panelHeight, y1));
                        y2 = Math.Max(0, Math.Min(panelHeight, y2));
                        g.DrawLine(pen, x1, y1, x2, y2);
                    }
                }

                // Plot water level
                using (Pen pen = new Pen(Color.Green, 2))
                {
                    float maxWaterLevel = (float)waterLevel.Max(Math.Abs);
                    if (maxWaterLevel == 0 || float.IsNaN(maxWaterLevel) || float.IsInfinity(maxWaterLevel)) maxWaterLevel = tidalAmplitude > 0 ? (float)tidalAmplitude : 1.0f;
                    for (int i = 0; i < gridPoints - 1; i++)
                    {
                        float x1 = i * panelWidth / (gridPoints - 1);
                        float x2 = (i + 1) * panelWidth / (gridPoints - 1);
                        float y1 = panelHeight / 2 - (float)(waterLevel[i] / maxWaterLevel) * panelHeight / 8;
                        float y2 = panelHeight / 2 - (float)(waterLevel[i + 1] / maxWaterLevel) * panelHeight / 8;
                        if (float.IsNaN(y1) || float.IsInfinity(y1) || float.IsNaN(y2) || float.IsInfinity(y2)) continue;
                        y1 = Math.Max(0, Math.Min(panelHeight, y1));
                        y2 = Math.Max(0, Math.Min(panelHeight, y2));
                        g.DrawLine(pen, x1, y1, x2, y2);
                    }
                }
            }

            else // 1D Simplified
            {
                double dx1D = estuaryLength / (gridPoints - 1);
                using (Pen pen = new Pen(Color.Blue, 2))
                {
                    float maxVelocity = (float)velocityProfile.Max(Math.Abs);
                    if (maxVelocity == 0 || float.IsNaN(maxVelocity) || float.IsInfinity(maxVelocity)) maxVelocity = 1.0f;
                    for (int i = 0; i < gridPoints - 1; i++)
                    {
                        float x1 = i * panelWidth / (gridPoints - 1);
                        float x2 = (i + 1) * panelWidth / (gridPoints - 1);
                        float y1 = panelHeight / 2 - (float)(velocityProfile[i] / maxVelocity) * panelHeight / 4;
                        float y2 = panelHeight / 2 - (float)(velocityProfile[i + 1] / maxVelocity) * panelHeight / 4;
                        if (float.IsNaN(y1) || float.IsInfinity(y1) || float.IsNaN(y2) || float.IsInfinity(y2)) continue;
                        y1 = Math.Max(0, Math.Min(panelHeight, y1));
                        y2 = Math.Max(0, Math.Min(panelHeight, y2));
                        g.DrawLine(pen, x1, y1, x2, y2);
                    }
                }

                using (Pen pen = new Pen(Color.Red, 2))
                {
                    float maxSalinity = (float)salinityProfile.Max();
                    if (maxSalinity == 0 || float.IsNaN(maxSalinity) || float.IsInfinity(maxSalinity)) maxSalinity = 35.0f;
                    for (int i = 0; i < gridPoints - 1; i++)
                    {
                        float x1 = i * panelWidth / (gridPoints - 1);
                        float x2 = (i + 1) * panelWidth / (gridPoints - 1);
                        float y1 = panelHeight - (float)(salinityProfile[i] / maxSalinity) * panelHeight / 4 - panelHeight / 2;
                        float y2 = panelHeight - (float)(salinityProfile[i + 1] / maxSalinity) * panelHeight / 4 - panelHeight / 2;
                        if (float.IsNaN(y1) || float.IsInfinity(y1) || float.IsNaN(y2) || float.IsInfinity(y2)) continue;
                        y1 = Math.Max(0, Math.Min(panelHeight, y1));
                        y2 = Math.Max(0, Math.Min(panelHeight, y2));
                        g.DrawLine(pen, x1, y1, x2, y2);
                    }
                }

                using (Pen pen = new Pen(Color.Green, 2))
                {
                    float maxWaterLevel = (float)waterLevel.Max(Math.Abs);
                    if (maxWaterLevel == 0 || float.IsNaN(maxWaterLevel) || float.IsInfinity(maxWaterLevel)) maxWaterLevel = tidalAmplitude > 0 ? (float)tidalAmplitude : 1.0f;
                    for (int i = 0; i < gridPoints - 1; i++)
                    {
                        float x1 = i * panelWidth / (gridPoints - 1);
                        float x2 = (i + 1) * panelWidth / (gridPoints - 1);
                        float y1 = panelHeight / 2 - (float)(waterLevel[i] / maxWaterLevel) * panelHeight / 8;
                        float y2 = panelHeight / 2 - (float)(waterLevel[i + 1] / maxWaterLevel) * panelHeight / 8;
                        if (float.IsNaN(y1) || float.IsInfinity(y1) || float.IsNaN(y2) || float.IsInfinity(y2)) continue;
                        y1 = Math.Max(0, Math.Min(panelHeight, y1));
                        y2 = Math.Max(0, Math.Min(panelHeight, y2));
                        g.DrawLine(pen, x1, y1, x2, y2);
                    }
                }
            }

            // Draw legend
            using (Font font = new Font("Verdana", 8.25F))
            {
                g.DrawString("Velocity (Blue)", font, Brushes.Blue, 10, 10);
                g.DrawString("Salinity (Red)", font, Brushes.Red, 10, 30);
                g.DrawString("Water Level (Green)", font, Brushes.Green, 10, 50);
                if (solverType == "2D Shallow Water" && enableWetDry)
                {
                    g.DrawString("Wet (Light Blue) / Dry (Tan)", font, Brushes.Black, 10, 70);
                }
            }
        }
    }
}