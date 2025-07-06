using System;
using System.Drawing;
using System.Windows.Forms;
using System.IO;
using MathNet.Numerics;
using MathNet.Numerics.IntegralTransforms;
using System.Numerics;

namespace EstuarineCirculationModeling
{
    public class LatticeBoltzmannLES : Form
    {
        private int gridSizeX = 200; // Lattice grid size (x)
        private int gridSizeY = 50;  // Lattice grid size (y)
        private double dx = 50.0;    // Spatial step (meters)
#pragma warning disable CS0414 // Suppress warning: dt is assigned but not used; reserved for future physical time scaling
        private readonly double dt = 0.1; // Time step (seconds)
#pragma warning restore CS0414
        private double[,,] density;   // Fluid density
        private double[,,] salinity;  // Salinity field (PSU)
        private double[,,] f;        // Distribution functions (D2Q9) for velocity
        private double[,,] fEq;      // Equilibrium distribution functions for velocity
        private double[,,] g;        // Distribution functions for salinity
        private double[,,] gEq;      // Equilibrium distribution functions for salinity
        private double[,] u;         // x-velocity (instantaneous)
        private double[,] v;         // y-velocity (instantaneous)
        private double[,] uAvg;      // Time-averaged x-velocity
        private double[,] vAvg;      // Time-averaged y-velocity
        private double[,] uSum;      // Sum for time-averaged velocity
        private double[,] vSum;      // Sum for time-averaged velocity
        private double[,] vorticity;  // Vorticity or Q-criterion field
        private double[,] smagorinskyConstant; // Smagorinsky constant for LES
        private double[,] localDx;   // Local spatial step for refinement
        private double[,] localTau;  // Local relaxation time for velocity
        private double[,] localTauSalinity; // Local relaxation time for salinity
        private readonly double cs2 = 1.0 / 3.0; // Speed of sound squared for D2Q9
        private readonly double[] weights = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 }; // D2Q9 weights
        private readonly int[,] directions = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 1 }, { -1, 1 }, { -1, -1 }, { 1, -1 } }; // D2Q9 directions
        private Timer simulationTimer;
        private bool isRunning;
        private double riverInflow = 0.1;     // m/s
        private double tidalAmplitude = 1.0;   // meters
        private double tidalPeriod = 43200.0;  // seconds
        private double smagorinskyC = 0.1;    // Smagorinsky constant
        private double reynoldsNumber = 1000.0; // Reynolds number
        private double oceanSalinity = 35.0;  // PSU
        private double riverSalinity = 0.0;   // PSU
        private double oceanTemperature = 20.0; // °C
        private string refinementMode = "None"; // Grid refinement mode
        private string filterWidthMode = "Lattice Spacing"; // Filter width mode
        private string vortexMode = "Vorticity Magnitude"; // Vorticity or Q-criterion
        private string vortexColorMode = "Velocity"; // Color by velocity or vorticity
        private TextBox gridSizeXTextBox;
        private TextBox gridSizeYTextBox;
        private TextBox dxTextBox;
        private TextBox riverInflowTextBox;
        private TextBox tidalAmplitudeTextBox;
        private TextBox tidalPeriodTextBox;
        private TextBox smagorinskyTextBox;
        private TextBox reynoldsNumberTextBox;
        private TextBox oceanSalinityTextBox;
        private TextBox riverSalinityTextBox;
        private TextBox oceanTemperatureTextBox;
        private ComboBox visualizationModeComboBox;
        private ComboBox refinementModeComboBox;
        private ComboBox filterWidthModeComboBox;
        private ComboBox vortexModeComboBox;
        private ComboBox vortexColorComboBox;
        private Button startButton;
        private Button pauseButton;
        private Button resetButton;
        private TextBox outputConsole;
        private Panel visualizationPanel;
        private Panel vorticityPanel;
        private string visualizationMode = "Velocity"; // Default visualization mode
        private int timeStepCount = 0; // For time averaging and spectrum output
        private int avgCount = 0;    // Number of timesteps for averaging

        public LatticeBoltzmannLES()
        {
            InitializeWindow();
            InitializeArrays();
            InitializeSimulation();
            InitializeTimer();
        }

        private void InitializeArrays()
        {
            density = new double[gridSizeX, gridSizeY, 1];
            salinity = new double[gridSizeX, gridSizeY, 1];
            f = new double[gridSizeX, gridSizeY, 9];
            fEq = new double[gridSizeX, gridSizeY, 9];
            g = new double[gridSizeX, gridSizeY, 9];
            gEq = new double[gridSizeX, gridSizeY, 9];
            u = new double[gridSizeX, gridSizeY];
            v = new double[gridSizeX, gridSizeY];
            uAvg = new double[gridSizeX, gridSizeY];
            vAvg = new double[gridSizeX, gridSizeY];
            uSum = new double[gridSizeX, gridSizeY];
            vSum = new double[gridSizeX, gridSizeY];
            vorticity = new double[gridSizeX, gridSizeY];
            smagorinskyConstant = new double[gridSizeX, gridSizeY];
            localDx = new double[gridSizeX, gridSizeY];
            localTau = new double[gridSizeX, gridSizeY];
            localTauSalinity = new double[gridSizeX, gridSizeY];
        }

        private void InitializeWindow()
        {
            this.Text = "Lattice-Boltzmann LES";
            this.ClientSize = new Size(1200, 600); // Increased width for dual panels
            this.FormBorderStyle = FormBorderStyle.FixedDialog;
            this.MaximizeBox = false;
            this.Font = new Font("Consolas", 9F);

            // Control Panel
            var controlPanel = new Panel
            {
                Location = new Point(12, 12),
                Size = new Size(300, 400),
                AutoScroll = true,
                BorderStyle = BorderStyle.FixedSingle
            };

            // Input Labels and TextBoxes
            int y = 20;
            var labels = new[] { "Grid Size X:", "Grid Size Y:", "Spatial Step (m):", "River Inflow (m/s):",
                                 "Tidal Amplitude (m):", "Tidal Period (s):", "Smagorinsky Constant:",
                                 "Reynolds Number:", "Ocean Salinity (PSU):", "River Salinity (PSU):",
                                 "Ocean Temperature (°C):" };
            var defaults = new[] { "200", "50", "50.0", "0.1", "1.0", "43200", "0.1", "1000.0", "35.0", "0.0", "20.0" };
            var textBoxes = new TextBox[11];
            for (int i = 0; i < 11; i++)
            {
                var label = new Label
                {
                    AutoSize = true,
                    Location = new Point(10, y),
                    Text = labels[i],
                    Font = new Font("Consolas", 9F)
                };
                y += 20;
                textBoxes[i] = new TextBox
                {
                    Location = new Point(10, y),
                    Size = new Size(250, 22),
                    Text = defaults[i],
                    BorderStyle = BorderStyle.FixedSingle,
                    Font = new Font("Consolas", 9F)
                };
                controlPanel.Controls.Add(label);
                controlPanel.Controls.Add(textBoxes[i]);
                y += 30;
            }
            gridSizeXTextBox = textBoxes[0];
            gridSizeYTextBox = textBoxes[1];
            dxTextBox = textBoxes[2];
            riverInflowTextBox = textBoxes[3];
            tidalAmplitudeTextBox = textBoxes[4];
            tidalPeriodTextBox = textBoxes[5];
            smagorinskyTextBox = textBoxes[6];
            reynoldsNumberTextBox = textBoxes[7];
            oceanSalinityTextBox = textBoxes[8];
            riverSalinityTextBox = textBoxes[9];
            oceanTemperatureTextBox = textBoxes[10];

            // Visualization Mode Dropdown
            visualizationModeComboBox = new ComboBox
            {
                Location = new Point(10, y),
                Size = new Size(250, 22),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = new Font("Consolas", 9F)
            };
            visualizationModeComboBox.Items.AddRange(new[] { "Velocity", "Salinity", "Time-Averaged Velocity" });
            visualizationModeComboBox.SelectedIndex = 0;
            visualizationModeComboBox.SelectedIndexChanged += (s, e) =>
            {
                visualizationMode = visualizationModeComboBox.SelectedItem.ToString();
                visualizationPanel.Invalidate();
            };
            controlPanel.Controls.Add(visualizationModeComboBox);
            y += 30;

            // Refinement Mode Dropdown
            refinementModeComboBox = new ComboBox
            {
                Location = new Point(10, y),
                Size = new Size(250, 22),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = new Font("Consolas", 9F)
            };
            refinementModeComboBox.Items.AddRange(new[] { "None", "River", "Tidal", "Both" });
            refinementModeComboBox.SelectedIndex = 0;
            refinementModeComboBox.SelectedIndexChanged += (s, e) =>
            {
                refinementMode = refinementModeComboBox.SelectedItem?.ToString() ?? "None";
                InitializeSimulation();
                visualizationPanel.Invalidate();
                vorticityPanel.Invalidate();
            };
            controlPanel.Controls.Add(refinementModeComboBox);
            y += 30;

            // Filter Width Mode Dropdown
            filterWidthModeComboBox = new ComboBox
            {
                Location = new Point(10, y),
                Size = new Size(250, 22),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = new Font("Consolas", 9F)
            };
            filterWidthModeComboBox.Items.AddRange(new[] { "Lattice Spacing", "Cell Area" });
            filterWidthModeComboBox.SelectedIndex = 0;
            filterWidthModeComboBox.SelectedIndexChanged += (s, e) =>
            {
                filterWidthMode = filterWidthModeComboBox.SelectedItem?.ToString() ?? "Lattice Spacing";
                InitializeSimulation();
                visualizationPanel.Invalidate();
                vorticityPanel.Invalidate();
            };
            controlPanel.Controls.Add(filterWidthModeComboBox);
            y += 30;

            // Vortex Mode Dropdown
            vortexModeComboBox = new ComboBox
            {
                Location = new Point(10, y),
                Size = new Size(250, 22),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = new Font("Consolas", 9F)
            };
            vortexModeComboBox.Items.AddRange(new[] { "Vorticity Magnitude", "Q-Criterion" });
            vortexModeComboBox.SelectedIndex = 0;
            vortexModeComboBox.SelectedIndexChanged += (s, e) =>
            {
                vortexMode = vortexModeComboBox.SelectedItem.ToString();
                vorticityPanel.Invalidate();
            };
            controlPanel.Controls.Add(vortexModeComboBox);
            y += 30;

            // Vortex Color Mode Dropdown
            vortexColorComboBox = new ComboBox
            {
                Location = new Point(10, y),
                Size = new Size(250, 22),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = new Font("Consolas", 9F)
            };
            vortexColorComboBox.Items.AddRange(new[] { "Velocity", "Vorticity" });
            vortexColorComboBox.SelectedIndex = 0;
            vortexColorComboBox.SelectedIndexChanged += (s, e) =>
            {
                vortexColorMode = vortexColorComboBox.SelectedItem.ToString();
                vorticityPanel.Invalidate();
            };
            controlPanel.Controls.Add(vortexColorComboBox);
            y += 30;

            // Buttons
            startButton = new Button
            {
                Location = new Point(10, y),
                Size = new Size(170, 25),
                Text = "Start",
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                Font = new Font("Consolas", 9F)
            };
            startButton.Click += StartButton_Click;
            y += 35;
            pauseButton = new Button
            {
                Location = new Point(10, y),
                Size = new Size(170, 25),
                Text = "Pause",
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                Enabled = false,
                Font = new Font("Consolas", 9F)
            };
            pauseButton.Click += PauseButton_Click;
            y += 35;
            resetButton = new Button
            {
                Location = new Point(10, y),
                Size = new Size(170, 25),
                Text = "Reset",
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black },
                Font = new Font("Consolas", 9F)
            };
            resetButton.Click += ResetButton_Click;
            controlPanel.Controls.Add(startButton);
            controlPanel.Controls.Add(pauseButton);
            controlPanel.Controls.Add(resetButton);

            // Visualization Panel (Original)
            visualizationPanel = new Panel
            {
                Location = new Point(318, 12),
                Size = new Size(400, 400),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            visualizationPanel.Paint += VisualizationPanel_Paint;

            // Vorticity Panel (New)
            vorticityPanel = new Panel
            {
                Location = new Point(724, 12),
                Size = new Size(400, 400),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            vorticityPanel.Paint += VorticityPanel_Paint;

            // Output Console
            outputConsole = new TextBox
            {
                Location = new Point(12, 420),
                Size = new Size(1176, 150),
                Multiline = true,
                ReadOnly = true,
                ScrollBars = ScrollBars.Vertical,
                BorderStyle = BorderStyle.FixedSingle,
                Font = new Font("Consolas", 9F)
            };

            this.Controls.Add(controlPanel);
            this.Controls.Add(visualizationPanel);
            this.Controls.Add(vorticityPanel);
            this.Controls.Add(outputConsole);
        }

        private void InitializeSimulation()
        {
            // Compute molecular viscosity from Reynolds number
            double characteristicLength = gridSizeX * dx;
            double nu = riverInflow * characteristicLength / reynoldsNumber;

            // Initialize localDx and localTau
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                {
                    double refinementFactor = 1.0;
                    if (refinementMode == "River" && i < gridSizeX / 10)
                        refinementFactor = 0.5; // 2x resolution near river
                    else if (refinementMode == "Tidal" && i > gridSizeX * 9 / 10)
                        refinementFactor = 0.5; // 2x resolution near tidal boundary
                    else if (refinementMode == "Both" && (i < gridSizeX / 10 || i > gridSizeX * 9 / 10))
                        refinementFactor = 0.5;
                    localDx[i, j] = dx * refinementFactor;

                    // Initialize localTau with molecular viscosity only
                    localTau[i, j] = 3.0 * nu * cs2 * (dx / localDx[i, j]) + 0.5;
                    localTauSalinity[i, j] = 3.0 * nu * cs2 * (dx / localDx[i, j]) + 0.5;
                }

            // Initialize density, salinity, and velocity fields
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                {
                    salinity[i, j, 0] = riverSalinity + (oceanSalinity - riverSalinity) * i / (gridSizeX - 1);
                    density[i, j, 0] = 1000.0 * (1.0 + 0.0008 * salinity[i, j, 0] - 0.00007 * (oceanTemperature - 20.0));
                    u[i, j] = riverInflow;
                    v[i, j] = 0.0;
                    uAvg[i, j] = 0.0;
                    vAvg[i, j] = 0.0;
                    uSum[i, j] = 0.0;
                    vSum[i, j] = 0.0;
                    vorticity[i, j] = 0.0;
                    smagorinskyConstant[i, j] = smagorinskyC;
                    for (int k = 0; k < 9; k++)
                    {
                        f[i, j, k] = weights[k] * density[i, j, 0];
                        fEq[i, j, k] = f[i, j, k];
                        g[i, j, k] = weights[k] * salinity[i, j, 0];
                        gEq[i, j, k] = g[i, j, k];
                    }
                }
            timeStepCount = 0;
            avgCount = 0;
        }

        private void InitializeTimer()
        {
            simulationTimer = new Timer
            {
                Interval = 100
            };
            simulationTimer.Tick += (s, e) => RunSimulationStep();
            isRunning = false;
            UpdateButtonStates();
        }

        private void StartButton_Click(object sender, EventArgs e)
        {
            try
            {
                gridSizeX = Math.Max(50, Math.Min(500, int.Parse(gridSizeXTextBox.Text)));
                gridSizeY = Math.Max(50, Math.Min(500, int.Parse(gridSizeYTextBox.Text)));
                dx = Math.Max(10.0, Math.Min(100.0, double.Parse(dxTextBox.Text)));
                riverInflow = Math.Max(0.01, Math.Min(10.0, double.Parse(riverInflowTextBox.Text)));
                tidalAmplitude = Math.Max(0.1, Math.Min(5.0, double.Parse(tidalAmplitudeTextBox.Text)));
                tidalPeriod = Math.Max(3600.0, Math.Min(86400.0, double.Parse(tidalPeriodTextBox.Text)));
                smagorinskyC = Math.Max(0.05, Math.Min(0.3, double.Parse(smagorinskyTextBox.Text)));
                reynoldsNumber = Math.Max(100.0, Math.Min(10000.0, double.Parse(reynoldsNumberTextBox.Text)));
                oceanSalinity = Math.Max(0.0, Math.Min(40.0, double.Parse(oceanSalinityTextBox.Text)));
                riverSalinity = Math.Max(0.0, Math.Min(40.0, double.Parse(riverSalinityTextBox.Text)));
                oceanTemperature = Math.Max(0.0, Math.Min(40.0, double.Parse(oceanTemperatureTextBox.Text)));
                refinementMode = refinementModeComboBox.SelectedItem?.ToString() ?? "None";
                filterWidthMode = filterWidthModeComboBox.SelectedItem?.ToString() ?? "Lattice Spacing";
                vortexMode = vortexModeComboBox.SelectedItem?.ToString() ?? "Vorticity Magnitude";
                vortexColorMode = vortexColorComboBox.SelectedItem?.ToString() ?? "Velocity";
                InitializeArrays();
                InitializeSimulation();
                simulationTimer.Start();
                isRunning = true;
                UpdateButtonStates();
                outputConsole.AppendText("LB-LES simulation started.\r\n");
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Error: {ex.Message}", "Input Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private void PauseButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isRunning = false;
            UpdateButtonStates();
            outputConsole.AppendText("LB-LES simulation paused.\r\n");
        }

        private void ResetButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isRunning = false;
            InitializeArrays();
            InitializeSimulation();
            visualizationPanel.Invalidate();
            vorticityPanel.Invalidate();
            outputConsole.Clear();
            outputConsole.AppendText("LB-LES simulation reset.\r\n");
            UpdateButtonStates();
        }

        private void UpdateButtonStates()
        {
            startButton.Enabled = !isRunning;
            pauseButton.Enabled = isRunning;
            resetButton.Enabled = true;
        }

        private void RunSimulationStep()
        {
            double currentTime = simulationTimer.Interval * 0.001 * Environment.TickCount / 1000.0;
            double tidalVelocity = tidalAmplitude * Math.Sin(2 * Math.PI * currentTime / tidalPeriod) / dx;
            timeStepCount++;
            avgCount++;

            // Streaming step for velocity
            double[,,] fTemp = new double[gridSizeX, gridSizeY, 9];
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                    for (int k = 0; k < 9; k++)
                    {
                        int nextI = (i + directions[k, 0] + gridSizeX) % gridSizeX;
                        int nextJ = (j + directions[k, 1] + gridSizeY) % gridSizeY;
                        fTemp[nextI, nextJ, k] = f[i, j, k];
                    }
            Array.Copy(fTemp, f, fTemp.Length);

            // Streaming step for salinity
            double[,,] gTemp = new double[gridSizeX, gridSizeY, 9];
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                    for (int k = 0; k < 9; k++)
                    {
                        int nextI = (i + directions[k, 0] + gridSizeX) % gridSizeX;
                        int nextJ = (j + directions[k, 1] + gridSizeY) % gridSizeY;
                        gTemp[nextI, nextJ, k] = g[i, j, k];
                    }
            Array.Copy(gTemp, g, gTemp.Length);

            // Compute molecular viscosity
            double characteristicLength = gridSizeX * dx;
            double nu = riverInflow * characteristicLength / reynoldsNumber;

            // Compute turbulence statistics
            double avgTKE = 0.0;
            double avgDissipation = 0.0;

            // Collision step with Smagorinsky LES and vorticity/Q-criterion calculation
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                {
                    // Compute macroscopic variables
                    density[i, j, 0] = 0.0;
                    u[i, j] = 0.0;
                    v[i, j] = 0.0;
                    salinity[i, j, 0] = 0.0;
                    for (int k = 0; k < 9; k++)
                    {
                        density[i, j, 0] += f[i, j, k];
                        u[i, j] += f[i, j, k] * directions[k, 0];
                        v[i, j] += f[i, j, k] * directions[k, 1];
                        salinity[i, j, 0] += g[i, j, k];
                    }
                    u[i, j] /= density[i, j, 0];
                    v[i, j] /= density[i, j, 0];
                    salinity[i, j, 0] /= density[i, j, 0];

                    // Update time-averaged velocities
                    uSum[i, j] += u[i, j];
                    vSum[i, j] += v[i, j];
                    uAvg[i, j] = uSum[i, j] / avgCount;
                    vAvg[i, j] = vSum[i, j] / avgCount;

                    // Update density based on salinity and temperature
                    density[i, j, 0] = 1000.0 * (1.0 + 0.0008 * salinity[i, j, 0] - 0.00007 * (oceanTemperature - 20.0));

                    // Compute strain rate and turbulent viscosity
                    double[,] S = new double[2, 2];
                    double S_magnitude = ComputeStrainRate(i, j, S);
                    double localCs = smagorinskyC * (1.0 - 0.5 * Math.Exp(-S_magnitude / 0.1));
                    double delta = filterWidthMode == "Lattice Spacing" ? localDx[i, j] : Math.Sqrt(localDx[i, j] * localDx[i, j]);
                    double nu_t = (localCs * delta) * (localCs * delta) * S_magnitude;
                    localTau[i, j] = 3.0 * (nu + nu_t) * cs2 * (dx / localDx[i, j]) + 0.5;
                    localTauSalinity[i, j] = 3.0 * (nu + nu_t) * cs2 * (dx / localDx[i, j]) + 0.5;

                    // Compute TKE
                    double uPrime = u[i, j] - uAvg[i, j];
                    double vPrime = v[i, j] - vAvg[i, j];
                    double tke = 0.5 * (uPrime * uPrime + vPrime * vPrime);
                    avgTKE += tke;

                    // Compute dissipation rate
                    double dissipation = 2.0 * nu_t * (S[0, 0] * S[0, 0] + S[0, 1] * S[0, 1] + S[1, 0] * S[1, 0] + S[1, 1] * S[1, 1]);
                    avgDissipation += dissipation;

                    // Compute vorticity or Q-criterion
                    vorticity[i, j] = ComputeVorticityOrQ(i, j, S);

                    // Compute equilibrium distributions for velocity
                    for (int k = 0; k < 9; k++)
                    {
                        double uDotC = directions[k, 0] * u[i, j] + directions[k, 1] * v[i, j];
                        double uSquared = u[i, j] * u[i, j] + v[i, j] * v[i, j];
                        fEq[i, j, k] = weights[k] * density[i, j, 0] * (1.0 + 3.0 * uDotC * (dx / localDx[i, j]) + 4.5 * uDotC * uDotC * (dx / localDx[i, j]) * (dx / localDx[i, j]) - 1.5 * uSquared * (dx / localDx[i, j]) * (dx / localDx[i, j]));
                        f[i, j, k] += (fEq[i, j, k] - f[i, j, k]) / localTau[i, j];
                    }

                    // Compute equilibrium distributions for salinity
                    for (int k = 0; k < 9; k++)
                    {
                        double uDotC = directions[k, 0] * u[i, j] + directions[k, 1] * v[i, j];
                        double uSquared = u[i, j] * u[i, j] + v[i, j] * v[i, j];
                        gEq[i, j, k] = weights[k] * salinity[i, j, 0] * density[i, j, 0] * (1.0 + 3.0 * uDotC * (dx / localDx[i, j]) + 4.5 * uDotC * uDotC * (dx / localDx[i, j]) * (dx / localDx[i, j]) - 1.5 * uSquared * (dx / localDx[i, j]) * (dx / localDx[i, j]));
                        g[i, j, k] += (gEq[i, j, k] - g[i, j, k]) / localTauSalinity[i, j];
                    }
                }

            // Average TKE and dissipation
            avgTKE /= (gridSizeX * gridSizeY);
            avgDissipation /= (gridSizeX * gridSizeY);

            // Compute energy spectrum every 100 timesteps
            if (timeStepCount % 100 == 0)
            {
                ComputeEnergySpectrum(currentTime);
            }

            // Apply boundary conditions
            for (int j = 0; j < gridSizeY; j++)
            {
                // River inflow (x=0)
                density[0, j, 0] = 1000.0 * (1.0 + 0.0008 * riverSalinity - 0.00007 * (oceanTemperature - 20.0));
                u[0, j] = riverInflow;
                v[0, j] = 0.0;
                salinity[0, j, 0] = riverSalinity;
                uSum[0, j] += u[0, j];
                vSum[0, j] += v[0, j];
                uAvg[0, j] = uSum[0, j] / avgCount;
                vAvg[0, j] = vSum[0, j] / avgCount;
                vorticity[0, j] = 0.0; // Simplified boundary condition
                for (int k = 0; k < 9; k++)
                {
                    f[0, j, k] = weights[k] * density[0, j, 0] * (1.0 + 3.0 * (directions[k, 0] * u[0, j] + directions[k, 1] * v[0, j]) * (dx / localDx[0, j]));
                    g[0, j, k] = weights[k] * salinity[0, j, 0] * density[0, j, 0];
                }

                // Tidal boundary (x=gridSizeX-1)
                density[gridSizeX - 1, j, 0] = 1000.0 * (1.0 + 0.0008 * oceanSalinity - 0.00007 * (oceanTemperature - 20.0));
                u[gridSizeX - 1, j] = tidalVelocity;
                v[gridSizeX - 1, j] = 0.0;
                salinity[gridSizeX - 1, j, 0] = oceanSalinity;
                uSum[gridSizeX - 1, j] += u[gridSizeX - 1, j];
                vSum[gridSizeX - 1, j] += v[gridSizeX - 1, j];
                uAvg[gridSizeX - 1, j] = uSum[gridSizeX - 1, j] / avgCount;
                vAvg[gridSizeX - 1, j] = vSum[gridSizeX - 1, j] / avgCount;
                vorticity[gridSizeX - 1, j] = 0.0; // Simplified boundary condition
                for (int k = 0; k < 9; k++)
                {
                    f[gridSizeX - 1, j, k] = weights[k] * density[gridSizeX - 1, j, 0] * (1.0 + 3.0 * (directions[k, 0] * u[gridSizeX - 1, j] + directions[k, 1] * v[gridSizeX - 1, j]) * (dx / localDx[gridSizeX - 1, j]));
                    g[gridSizeX - 1, j, k] = weights[k] * salinity[gridSizeX - 1, j, 0] * density[gridSizeX - 1, j, 0];
                }
            }

            visualizationPanel.Invalidate();
            vorticityPanel.Invalidate();
            UpdateOutputConsole(currentTime, avgTKE, avgDissipation);
        }

        private double ComputeStrainRate(int i, int j, double[,] S)
        {
            for (int k = 0; k < 9; k++)
            {
                double fNeq = f[i, j, k] - fEq[i, j, k];
                double cx = directions[k, 0] * (dx / localDx[i, j]);
                double cy = directions[k, 1] * (dx / localDx[i, j]);
                S[0, 0] += fNeq * cx * cx;
                S[0, 1] += fNeq * cx * cy;
                S[1, 0] += fNeq * cy * cx;
                S[1, 1] += fNeq * cy * cy;
            }
            double factor = -1.0 / (2.0 * density[i, j, 0] * cs2 * localTau[i, j]);
            S[0, 0] *= factor;
            S[0, 1] *= factor;
            S[1, 0] *= factor;
            S[1, 1] *= factor;
            return Math.Sqrt(2.0 * (S[0, 0] * S[0, 0] + S[0, 1] * S[0, 1] + S[1, 0] * S[1, 0] + S[1, 1] * S[1, 1]));
        }

        private double ComputeVorticityOrQ(int i, int j, double[,] S)
        {
            // Compute velocity gradients using finite differences
            double dudx, dudy, dvdx, dvdy;
            if (i == 0 || i == gridSizeX - 1 || j == 0 || j == gridSizeY - 1)
            {
                dudx = dudy = dvdx = dvdy = 0.0; // Simplified at boundaries
            }
            else
            {
                dudx = (u[i + 1, j] - u[i - 1, j]) / (2.0 * localDx[i, j]);
                dudy = (u[i, j + 1] - u[i, j - 1]) / (2.0 * localDx[i, j]);
                dvdx = (v[i + 1, j] - v[i - 1, j]) / (2.0 * localDx[i, j]);
                dvdy = (v[i, j + 1] - v[i, j - 1]) / (2.0 * localDx[i, j]);
            }

            if (vortexMode == "Vorticity Magnitude")
            {
                // Vorticity: |dv/dx - du/dy|
                return Math.Abs(dvdx - dudy);
            }
            else // Q-Criterion
            {
                // Rotation rate tensor Omega
                double[,] Omega = new double[2, 2];
                Omega[0, 1] = 0.5 * (dudx - dvdy); // Omega_xy = 0.5 * (du/dx - dv/dy)
                Omega[1, 0] = -Omega[0, 1]; // Omega_yx = -Omega_xy
                Omega[0, 0] = 0.0;
                Omega[1, 1] = 0.0;

                // Q = 0.5 * (||Omega||^2 - ||S||^2)
                double OmegaNorm = Omega[0, 1] * Omega[0, 1] + Omega[1, 0] * Omega[1, 0];
                double SNorm = S[0, 0] * S[0, 0] + S[0, 1] * S[0, 1] + S[1, 0] * S[1, 0] + S[1, 1] * S[1, 1];
                return 0.5 * (OmegaNorm - SNorm);
            }
        }

        private void ComputeEnergySpectrum(double currentTime)
        {
            // Extract velocity at mid-y plane
            int midY = gridSizeY / 2;
            Complex[] uData = new Complex[gridSizeX];
            Complex[] vData = new Complex[gridSizeX];
            for (int i = 0; i < gridSizeX; i++)
            {
                uData[i] = new Complex(u[i, midY], 0.0);
                vData[i] = new Complex(v[i, midY], 0.0);
            }

            // Compute FFT
            Fourier.Forward(uData, FourierOptions.Matlab);
            Fourier.Forward(vData, FourierOptions.Matlab);

            // Compute power spectral density
            double[] Ek = new double[gridSizeX / 2];
            double dk = 2.0 * Math.PI / (gridSizeX * dx); // Wavenumber increment
            for (int i = 0; i < gridSizeX / 2; i++)
            {
                Ek[i] = 0.5 * (uData[i].MagnitudeSquared() + vData[i].MagnitudeSquared());
            }

            // Output to file
            try
            {
                using (StreamWriter writer = new StreamWriter("energy_spectrum.txt", true))
                {
                    writer.WriteLine($"Time: {currentTime:F2}s");
                    writer.WriteLine("k\tE(k)");
                    for (int i = 0; i < gridSizeX / 2; i++)
                    {
                        double k = i * dk;
                        writer.WriteLine($"{k:F6}\t{Ek[i]:F6}");
                    }
                    writer.WriteLine();
                }
            }
            catch (Exception ex)
            {
                outputConsole.AppendText($"Error writing energy spectrum: {ex.Message}\r\n");
            }
        }

        private void VisualizationPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            float cellWidth = (float)visualizationPanel.Width / gridSizeX;
            float cellHeight = (float)visualizationPanel.Height / gridSizeY;

            if (visualizationMode == "Velocity")
            {
                double maxVelocity = 0.0;
                for (int i = 0; i < gridSizeX; i++)
                    for (int j = 0; j < gridSizeY; j++)
                    {
                        double velocity = Math.Sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j]);
                        maxVelocity = Math.Max(maxVelocity, velocity);
                    }
                maxVelocity = Math.Max(maxVelocity, 0.1);

                for (int i = 0; i < gridSizeX; i++)
                    for (int j = 0; j < gridSizeY; j++)
                    {
                        double velocity = Math.Sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j]);
                        int colorValue = (int)(255 * velocity / maxVelocity);
                        using (Brush brush = new SolidBrush(Color.FromArgb(colorValue, 0, 255 - colorValue)))
                        {
                            g.FillRectangle(brush, i * cellWidth, j * cellHeight, cellWidth, cellHeight);
                        }
                        if ((refinementMode == "River" && i < gridSizeX / 10) ||
                            (refinementMode == "Tidal" && i > gridSizeX * 9 / 10) ||
                            (refinementMode == "Both" && (i < gridSizeX / 10 || i > gridSizeX * 9 / 10)))
                        {
                            g.DrawRectangle(Pens.Yellow, i * cellWidth, j * cellHeight, cellWidth, cellHeight);
                        }
                    }

                for (int y = 0; y < 100; y++)
                {
                    int colorValue = (int)(255 * y / 100.0);
                    using (Brush brush = new SolidBrush(Color.FromArgb(colorValue, 0, 255 - colorValue)))
                    {
                        g.FillRectangle(brush, visualizationPanel.Width - 20, y * 4, 10, 4);
                    }
                }
                g.DrawString("Velocity", new Font("Consolas", 8), Brushes.Black, visualizationPanel.Width - 30, 100);
            }
            else if (visualizationMode == "Salinity")
            {
                double minSalinity = double.MaxValue, maxSalinity = double.MinValue;
                for (int i = 0; i < gridSizeX; i++)
                    for (int j = 0; j < gridSizeY; j++)
                    {
                        minSalinity = Math.Min(minSalinity, salinity[i, j, 0]);
                        maxSalinity = Math.Max(maxSalinity, salinity[i, j, 0]);
                    }
                maxSalinity = Math.Max(maxSalinity, minSalinity + 0.1);

                for (int i = 0; i < gridSizeX; i++)
                    for (int j = 0; j < gridSizeY; j++)
                    {
                        int colorValue = (int)(255 * (salinity[i, j, 0] - minSalinity) / (maxSalinity - minSalinity));
                        using (Brush brush = new SolidBrush(Color.FromArgb(0, colorValue, 255 - colorValue)))
                        {
                            g.FillRectangle(brush, i * cellWidth, j * cellHeight, cellWidth, cellHeight);
                        }
                        if ((refinementMode == "River" && i < gridSizeX / 10) ||
                            (refinementMode == "Tidal" && i > gridSizeX * 9 / 10) ||
                            (refinementMode == "Both" && (i < gridSizeX / 10 || i > gridSizeX * 9 / 10)))
                        {
                            g.DrawRectangle(Pens.Yellow, i * cellWidth, j * cellHeight, cellWidth, cellHeight);
                        }
                    }

                for (int y = 0; y < 100; y++)
                {
                    int colorValue = (int)(255 * y / 100.0);
                    using (Brush brush = new SolidBrush(Color.FromArgb(0, colorValue, 255 - colorValue)))
                    {
                        g.FillRectangle(brush, visualizationPanel.Width - 20, y * 4, 10, 4);
                    }
                }
                g.DrawString("Salinity", new Font("Consolas", 8), Brushes.Black, visualizationPanel.Width - 30, 100);
            }
            else // Time-Averaged Velocity
            {
                double maxVelocity = 0.0;
                for (int i = 0; i < gridSizeX; i++)
                    for (int j = 0; j < gridSizeY; j++)
                    {
                        double velocity = Math.Sqrt(uAvg[i, j] * uAvg[i, j] + vAvg[i, j] * vAvg[i, j]);
                        maxVelocity = Math.Max(maxVelocity, velocity);
                    }
                maxVelocity = Math.Max(maxVelocity, 0.1);

                for (int i = 0; i < gridSizeX; i++)
                    for (int j = 0; j < gridSizeY; j++)
                    {
                        double velocity = Math.Sqrt(uAvg[i, j] * uAvg[i, j] + vAvg[i, j] * vAvg[i, j]);
                        int colorValue = (int)(255 * velocity / maxVelocity);
                        using (Brush brush = new SolidBrush(Color.FromArgb(colorValue, 0, 255 - colorValue)))
                        {
                            g.FillRectangle(brush, i * cellWidth, j * cellHeight, cellWidth, cellHeight);
                        }
                        if ((refinementMode == "River" && i < gridSizeX / 10) ||
                            (refinementMode == "Tidal" && i > gridSizeX * 9 / 10) ||
                            (refinementMode == "Both" && (i < gridSizeX / 10 || i > gridSizeX * 9 / 10)))
                        {
                            g.DrawRectangle(Pens.Yellow, i * cellWidth, j * cellHeight, cellWidth, cellHeight);
                        }
                    }

                for (int y = 0; y < 100; y++)
                {
                    int colorValue = (int)(255 * y / 100.0);
                    using (Brush brush = new SolidBrush(Color.FromArgb(colorValue, 0, 255 - colorValue)))
                    {
                        g.FillRectangle(brush, visualizationPanel.Width - 20, y * 4, 10, 4);
                    }
                }
                g.DrawString("Time-Avg Velocity", new Font("Consolas", 8), Brushes.Black, visualizationPanel.Width - 30, 100);
            }
        }

        private void VorticityPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            float cellWidth = (float)vorticityPanel.Width / gridSizeX;
            float cellHeight = (float)vorticityPanel.Height / gridSizeY;

            // Find max vorticity/Q-criterion and max color value (velocity or vorticity)
            double maxVorticity = 0.0, maxColorValue = 0.0;
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                {
                    maxVorticity = Math.Max(maxVorticity, Math.Abs(vorticity[i, j]));
                    if (vortexColorMode == "Velocity")
                    {
                        maxColorValue = Math.Max(maxColorValue, Math.Sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j]));
                    }
                    else
                    {
                        maxColorValue = Math.Max(maxColorValue, Math.Abs(vorticity[i, j]));
                    }
                }
            maxVorticity = Math.Max(maxVorticity, 0.01);
            maxColorValue = Math.Max(maxColorValue, 0.1);

            // Threshold for isosurfaces (10% of max)
            double threshold = vortexMode == "Q-Criterion" ? 0.0 : 0.1 * maxVorticity;

            // Draw isosurfaces
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                {
                    if ((vortexMode == "Q-Criterion" && vorticity[i, j] > threshold) ||
                        (vortexMode == "Vorticity Magnitude" && vorticity[i, j] > threshold))
                    {
                        double colorScalar = vortexColorMode == "Velocity"
                            ? Math.Sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j])
                            : vorticity[i, j];
                        int colorValue = (int)(255 * colorScalar / maxColorValue);
                        using (Brush brush = new SolidBrush(Color.FromArgb(colorValue, 0, 255 - colorValue)))
                        {
                            g.FillEllipse(brush, i * cellWidth, j * cellHeight, cellWidth, cellHeight);
                        }
                    }
                    if ((refinementMode == "River" && i < gridSizeX / 10) ||
                        (refinementMode == "Tidal" && i > gridSizeX * 9 / 10) ||
                        (refinementMode == "Both" && (i < gridSizeX / 10 || i > gridSizeX * 9 / 10)))
                    {
                        g.DrawRectangle(Pens.Yellow, i * cellWidth, j * cellHeight, cellWidth, cellHeight);
                    }
                }

            // Color bar
            for (int y = 0; y < 100; y++)
            {
                int colorValue = (int)(255 * y / 100.0);
                using (Brush brush = new SolidBrush(Color.FromArgb(colorValue, 0, 255 - colorValue)))
                {
                    g.FillRectangle(brush, vorticityPanel.Width - 20, y * 4, 10, 4);
                }
            }
            g.DrawString(vortexColorMode, new Font("Consolas", 8), Brushes.Black, vorticityPanel.Width - 30, 100);
        }

        private void UpdateOutputConsole(double time, double avgTKE, double avgDissipation)
        {
            double avgVelocity = 0.0, avgSalinity = 0.0, avgVorticity = 0.0;
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                {
                    avgVelocity += Math.Sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j]);
                    avgSalinity += salinity[i, j, 0];
                    avgVorticity += Math.Abs(vorticity[i, j]);
                }
            avgVelocity /= (gridSizeX * gridSizeY);
            avgSalinity /= (gridSizeX * gridSizeY);
            avgVorticity /= (gridSizeX * gridSizeY);
            outputConsole.AppendText($"Time: {time:F2}s | Avg Velocity: {avgVelocity:F4} m/s | Avg Salinity: {avgSalinity:F2} PSU | Avg TKE: {avgTKE:F6} m²/s² | Avg Dissipation: {avgDissipation:F6} m²/s³ | Avg Vorticity: {avgVorticity:F6}\r\n");
        }

        public static void ShowLBLESWindow()
        {
            using (var form = new LatticeBoltzmannLES())
            {
                form.ShowDialog();
            }
        }

        protected override void Dispose(bool disposing)
        {
            if (disposing)
            {
                simulationTimer?.Dispose();
            }
            base.Dispose(disposing);
        }
    }
}