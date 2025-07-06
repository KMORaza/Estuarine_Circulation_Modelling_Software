using System;
using System.Drawing;
using System.Windows.Forms;

namespace EstuarineCirculationModeling
{
    public class LatticeBoltzmannLES : Form
    {
        private readonly int gridSizeX = 200; // Lattice grid size (x)
        private readonly int gridSizeY = 50;  // Lattice grid size (y)
        private readonly double dx = 50.0;    // Spatial step (meters)
        private readonly double dt = 0.1;     // Time step (seconds)
        private readonly double[,,] density;   // Fluid density
        private readonly double[,,] f;        // Distribution functions (D2Q9)
        private readonly double[,,] fEq;      // Equilibrium distribution functions
        private readonly double[,] u;         // x-velocity
        private readonly double[,] v;         // y-velocity
        private readonly double[,] smagorinskyConstant; // Smagorinsky constant for LES
        private readonly double tau = 0.8;    // Relaxation time
        private readonly double cs2 = 1.0 / 3.0; // Speed of sound squared for D2Q9
        private readonly double[] weights = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 }; // D2Q9 weights
        private readonly int[,] directions = { { 0, 0 }, { 1, 0 }, { 0, 1 }, { -1, 0 }, { 0, -1 }, { 1, 1 }, { -1, 1 }, { -1, -1 }, { 1, -1 } }; // D2Q9 directions
        private Timer simulationTimer;
        private bool isRunning;
        private double riverInflow = 0.1;     // m/s
        private double tidalAmplitude = 1.0;   // meters
        private double tidalPeriod = 43200.0;  // seconds
        private double smagorinskyC = 0.1;    // Smagorinsky constant
        private Panel visualizationPanel;
        private TextBox riverInflowTextBox;
        private TextBox tidalAmplitudeTextBox;
        private TextBox tidalPeriodTextBox;
        private TextBox smagorinskyTextBox;
        private Button startButton;
        private Button pauseButton;
        private Button resetButton;
        private TextBox outputConsole;

        public LatticeBoltzmannLES()
        {
            density = new double[gridSizeX, gridSizeY, 1]; // 3D array for scalar density
            f = new double[gridSizeX, gridSizeY, 9];      // 3D array for D2Q9
            fEq = new double[gridSizeX, gridSizeY, 9];     // 3D array for equilibrium
            u = new double[gridSizeX, gridSizeY];
            v = new double[gridSizeX, gridSizeY];
            smagorinskyConstant = new double[gridSizeX, gridSizeY];
            InitializeWindow();
            InitializeSimulation();
            InitializeTimer();
        }

        private void InitializeWindow()
        {
            this.Text = "Lattice-Boltzmann LES";
            this.ClientSize = new Size(900, 600);
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
            var labels = new[] { "River Inflow (m/s):", "Tidal Amplitude (m):", "Tidal Period (s):", "Smagorinsky Constant:" };
            var defaults = new[] { "0.1", "1.0", "43200", "0.1" };
            var textBoxes = new TextBox[4];
            for (int i = 0; i < 4; i++)
            {
                var label = new Label
                {
                    AutoSize = true,
                    Location = new Point(10, y),
                    Text = labels[i]
                };
                y += 20;
                textBoxes[i] = new TextBox
                {
                    Location = new Point(10, y),
                    Size = new Size(250, 22),
                    Text = defaults[i],
                    BorderStyle = BorderStyle.FixedSingle
                };
                controlPanel.Controls.Add(label);
                controlPanel.Controls.Add(textBoxes[i]);
                y += 30;
            }
            riverInflowTextBox = textBoxes[0];
            tidalAmplitudeTextBox = textBoxes[1];
            tidalPeriodTextBox = textBoxes[2];
            smagorinskyTextBox = textBoxes[3];

            // Buttons
            startButton = new Button
            {
                Location = new Point(10, y),
                Size = new Size(170, 25),
                Text = "Start",
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black }
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
                Enabled = false
            };
            pauseButton.Click += PauseButton_Click;
            y += 35;
            resetButton = new Button
            {
                Location = new Point(10, y),
                Size = new Size(170, 25),
                Text = "Reset",
                FlatStyle = FlatStyle.Flat,
                FlatAppearance = { BorderSize = 1, BorderColor = Color.Black }
            };
            resetButton.Click += ResetButton_Click;
            controlPanel.Controls.Add(startButton);
            controlPanel.Controls.Add(pauseButton);
            controlPanel.Controls.Add(resetButton);

            // Visualization Panel
            visualizationPanel = new Panel
            {
                Location = new Point(318, 12),
                Size = new Size(560, 400),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            visualizationPanel.Paint += VisualizationPanel_Paint;

            outputConsole = new TextBox
            {
                Location = new Point(12, 420),
                Size = new Size(860, 150),
                Multiline = true,
                ReadOnly = true,
                ScrollBars = ScrollBars.Vertical,
                BorderStyle = BorderStyle.FixedSingle
            };

            // LatticeBoltxmannLES
            this.Controls.Add(controlPanel);
            this.Controls.Add(visualizationPanel);
            this.Controls.Add(outputConsole);
        }

        private void InitializeSimulation()
        {
            // Initialize density and velocity fields
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                {
                    density[i, j, 0] = 1000.0; // kg/m³
                    u[i, j] = riverInflow;
                    v[i, j] = 0.0;
                    smagorinskyConstant[i, j] = smagorinskyC;
                    for (int k = 0; k < 9; k++)
                    {
                        f[i, j, k] = weights[k] * density[i, j, 0];
                        fEq[i, j, k] = f[i, j, k];
                    }
                }
        }

        private void InitializeTimer()
        {
            simulationTimer = new Timer
            {
                Interval = 100 // Update every 100ms
            };
            simulationTimer.Tick += (s, e) => RunSimulationStep();
            isRunning = false;
            UpdateButtonStates();
        }

        private void StartButton_Click(object sender, EventArgs e)
        {
            try
            {
                riverInflow = Math.Max(0.01, Math.Min(10.0, double.Parse(riverInflowTextBox.Text)));
                tidalAmplitude = Math.Max(0.1, Math.Min(5.0, double.Parse(tidalAmplitudeTextBox.Text)));
                tidalPeriod = Math.Max(3600.0, Math.Min(86400.0, double.Parse(tidalPeriodTextBox.Text)));
                smagorinskyC = Math.Max(0.05, Math.Min(0.3, double.Parse(smagorinskyTextBox.Text)));
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
            InitializeSimulation();
            visualizationPanel.Invalidate();
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
            // Apply tidal boundary condition at x = gridSizeX - 1
            double currentTime = simulationTimer.Interval * 0.001 * Environment.TickCount / 1000.0;
            double tidalVelocity = tidalAmplitude * Math.Sin(2 * Math.PI * currentTime / tidalPeriod) / dx;

            // Streaming step
            double[,,] fTemp = new double[gridSizeX, gridSizeY, 9];
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                    for (int k = 0; k < 9; k++)
                    {
                        int nextI = (i + directions[k, 0] + gridSizeX) % gridSizeX;
                        int nextJ = (j + directions[k, 1] + gridSizeY) % gridSizeY;
                        fTemp[nextI, nextJ, k] = f[i, j, k];
                    }

            // Update distribution functions
            Array.Copy(fTemp, f, fTemp.Length);

            // Collision step with Smagorinsky LES
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                {
                    // Compute macroscopic variables
                    density[i, j, 0] = 0.0;
                    u[i, j] = 0.0;
                    v[i, j] = 0.0;
                    for (int k = 0; k < 9; k++)
                    {
                        density[i, j, 0] += f[i, j, k];
                        u[i, j] += f[i, j, k] * directions[k, 0];
                        v[i, j] += f[i, j, k] * directions[k, 1];
                    }
                    u[i, j] /= density[i, j, 0];
                    v[i, j] /= density[i, j, 0];

                    // Compute strain rate tensor for Smagorinsky model
                    double S = ComputeStrainRate(i, j);
                    double tauTotal = tau + 3.0 * smagorinskyC * smagorinskyC * S * dx * dx / dt;

                    // Compute equilibrium distribution
                    for (int k = 0; k < 9; k++)
                    {
                        double uDotC = directions[k, 0] * u[i, j] + directions[k, 1] * v[i, j];
                        double uSquared = u[i, j] * u[i, j] + v[i, j] * v[i, j];
                        fEq[i, j, k] = weights[k] * density[i, j, 0] * (1.0 + 3.0 * uDotC + 4.5 * uDotC * uDotC - 1.5 * uSquared);
                        f[i, j, k] += (fEq[i, j, k] - f[i, j, k]) / tauTotal;
                    }
                }

            // Apply boundary conditions (inflow at x=0, tidal at x=gridSizeX-1)
            for (int j = 0; j < gridSizeY; j++)
            {
                // River inflow
                density[0, j, 0] = 1000.0;
                u[0, j] = riverInflow;
                v[0, j] = 0.0;
                for (int k = 0; k < 9; k++)
                    f[0, j, k] = weights[k] * density[0, j, 0] * (1.0 + 3.0 * (directions[k, 0] * u[0, j] + directions[k, 1] * v[0, j]));

                // Tidal boundary
                density[gridSizeX - 1, j, 0] = 1000.0;
                u[gridSizeX - 1, j] = tidalVelocity;
                v[gridSizeX - 1, j] = 0.0;
                for (int k = 0; k < 9; k++)
                    f[gridSizeX - 1, j, k] = weights[k] * density[gridSizeX - 1, j, 0] * (1.0 + 3.0 * (directions[k, 0] * u[gridSizeX - 1, j] + directions[k, 1] * v[gridSizeX - 1, j]));
            }

            visualizationPanel.Invalidate();
            UpdateOutputConsole(currentTime);
        }

        private double ComputeStrainRate(int i, int j)
        {
            // Compute strain rate tensor S_{\alpha\beta} using non-equilibrium moments
            double[,] S = new double[2, 2]; // 2x2 strain rate tensor for 2D
            for (int k = 0; k < 9; k++)
            {
                double fNeq = f[i, j, k] - fEq[i, j, k]; // Non-equilibrium part
                double cx = directions[k, 0];
                double cy = directions[k, 1];
                S[0, 0] += fNeq * cx * cx; // S_xx
                S[0, 1] += fNeq * cx * cy; // S_xy
                S[1, 0] += fNeq * cy * cx; // S_yx (symmetric to S_xy)
                S[1, 1] += fNeq * cy * cy; // S_yy
            }

            // Apply scaling factor: S_{\alpha\beta} = -1/(2 \rho c_s^2 \tau) * \sum f_i^{neq} c_{i\alpha} c_{i\beta}
            double factor = -1.0 / (2.0 * density[i, j, 0] * cs2 * tau);
            S[0, 0] *= factor;
            S[0, 1] *= factor;
            S[1, 0] *= factor;
            S[1, 1] *= factor;

            // Compute strain rate magnitude |S| = sqrt(2 S_{\alpha\beta} S_{\alpha\beta})
            double magnitude = Math.Sqrt(2.0 * (S[0, 0] * S[0, 0] + S[0, 1] * S[0, 1] + S[1, 0] * S[1, 0] + S[1, 1] * S[1, 1]));
            return magnitude;
        }

        private void VisualizationPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            float cellWidth = (float)visualizationPanel.Width / gridSizeX;
            float cellHeight = (float)visualizationPanel.Height / gridSizeY;

            // Compute max velocity for color scaling
            double maxVelocity = 0.0;
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                {
                    double velocity = Math.Sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j]);
                    maxVelocity = Math.Max(maxVelocity, velocity);
                }
            maxVelocity = Math.Max(maxVelocity, 0.1); // Avoid division by zero

            // Draw velocity field
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                {
                    double velocity = Math.Sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j]);
                    int colorValue = (int)(255 * velocity / maxVelocity);
                    using (Brush brush = new SolidBrush(Color.FromArgb(colorValue, 0, 255 - colorValue)))
                    {
                        g.FillRectangle(brush, i * cellWidth, j * cellHeight, cellWidth, cellHeight);
                    }
                }
        }

        private void UpdateOutputConsole(double time)
        {
            double avgVelocity = 0.0;
            for (int i = 0; i < gridSizeX; i++)
                for (int j = 0; j < gridSizeY; j++)
                    avgVelocity += Math.Sqrt(u[i, j] * u[i, j] + v[i, j] * v[i, j]);
            avgVelocity /= (gridSizeX * gridSizeY);
            outputConsole.AppendText($"Time: {time:F2}s | Avg Velocity: {avgVelocity:F4} m/s\r\n");
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