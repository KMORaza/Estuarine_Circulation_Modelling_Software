using System;
using System.Drawing;
using System.Windows.Forms;

namespace EstuarineCirculationModeling
{
    public class TurbulenceAndMixing
    {
        private Form turbulenceWindow;
        private PictureBox visualizationBox;
        private TextBox outputTextBox;
        private Button startButton;
        private Button pauseButton;
        private Button resetButton;
        private TextBox mixingCoefficientTextBox;
        private TextBox tidalAmplitudeTextBox;
        private TextBox estuaryDepthTextBox;
        private ComboBox sliceSelector;
        private Timer simulationTimer;
        private double[,] salinity; // 2D salinity field (x, z)
        private double[,] u; // Horizontal velocity (m/s)
        private double[,] w; // Vertical velocity (m/s)
        private double[,] p_h; // Hydrostatic pressure (Pa)
        private double[,] p_nh; // Non-hydrostatic pressure (Pa)
        private double[,] k; // Turbulent kinetic energy
        private double[,] epsilon; // Dissipation rate
        private double[,] buoyancy; // Buoyancy field (m/s²)
        private double[,] K_x; // Horizontal turbulent mixing coefficient (m²/s)
        private double[,] K_z; // Vertical turbulent mixing coefficient (m²/s)
        private double mixingCoefficient; // Background mixing coefficient
        private double tidalAmplitude; // Tidal amplitude (m)
        private double estuaryDepth; // Estuary depth (m)
        private double estuaryLength = 10000.0; // Fixed estuary length (m)
        private int gridPointsX = 100; // Horizontal grid points
        private int gridPointsZ = 50; // Vertical grid points
        private double dx; // Horizontal grid spacing
        private double dz; // Vertical grid spacing
        private double time; // Simulation time
        private double dt; // Adaptive time step
        private bool isRunning;
        private string currentSlice = "Surface"; // Default visualization
        private const double C_mu = 0.09; // k-ε model constant
        private const double sigma_k = 1.0; // k-ε model constant
        private const double sigma_epsilon = 1.3; // k-ε model constant
        private const double C1_epsilon = 1.44; // k-ε model constant
        private const double C2_epsilon = 1.92; // k-ε model constant
        private const double g = 9.81; // Gravitational acceleration
        private const double rho_0 = 1000.0; // Reference density (kg/m³)
        private const double beta = 0.8; // Saline contraction coefficient (kg/m³ per PSU)
        private const double refSalinity = 35.0; // Reference salinity (PSU)
        private const double minK = 1e-6; // Minimum TKE
        private const double minEpsilon = 1e-8; // Minimum dissipation rate
        private const double maxSalinity = 35.0; // Maximum salinity (PSU)
        private const double vonKarman = 0.41; // von Karman constant for logarithmic profile
        private const double z0 = 0.001; // Roughness length (m)
        private const double pressureRelaxation = 0.5; // Relaxation factor for pressure solver
        private const double minBuoyancy = -0.3; // Min buoyancy for visualization (m/s²)
        private const double maxBuoyancy = 0.0; // Max buoyancy for visualization (m/s²)
        private const double cflFactor = 0.7; // CFL factor for time step
        private const double minMixing = 1e-4; // Minimum mixing coefficient (m²/s)
        private const double maxMixing = 0.1; // Maximum mixing coefficient (m²/s)
        private readonly Font verdanaFont = new Font("Verdana", 8.25F);
        private readonly Random random = new Random(); // For initial salinity perturbation

        public TurbulenceAndMixing()
        {
            InitializeWindow();
            InitializeSimulation();
        }

        private void InitializeWindow()
        {
            turbulenceWindow = new Form
            {
                Text = "Turbulence and Mixing",
                Size = new Size(800, 600),
                FormBorderStyle = FormBorderStyle.FixedDialog,
                MaximizeBox = false,
                Font = verdanaFont
            };

            // Control panel
            Panel controlPanel = new Panel
            {
                Location = new Point(10, 10),
                Size = new Size(250, 540),
                BorderStyle = BorderStyle.FixedSingle,
                AutoScroll = true,
                Font = verdanaFont
            };

            // Labels and TextBoxes
            Label mixingLabel = new Label
            {
                Location = new Point(10, 10),
                Size = new Size(200, 20),
                Text = "Background Mixing (m²/s):",
                Font = verdanaFont
            };
            mixingCoefficientTextBox = new TextBox
            {
                Location = new Point(10, 30),
                Size = new Size(200, 20),
                Text = "0.01",
                Font = verdanaFont
            };
            Label tidalLabel = new Label
            {
                Location = new Point(10, 60),
                Size = new Size(200, 20),
                Text = "Tidal Amplitude (m):",
                Font = verdanaFont
            };
            tidalAmplitudeTextBox = new TextBox
            {
                Location = new Point(10, 80),
                Size = new Size(200, 20),
                Text = "2.0",
                Font = verdanaFont
            };
            Label depthLabel = new Label
            {
                Location = new Point(10, 110),
                Size = new Size(200, 20),
                Text = "Estuary Depth (m):",
                Font = verdanaFont
            };
            estuaryDepthTextBox = new TextBox
            {
                Location = new Point(10, 130),
                Size = new Size(200, 20),
                Text = "10.0",
                Font = verdanaFont
            };
            Label sliceLabel = new Label
            {
                Location = new Point(10, 160),
                Size = new Size(200, 20),
                Text = "Visualization Slice:",
                Font = verdanaFont
            };
            sliceSelector = new ComboBox
            {
                Location = new Point(10, 180),
                Size = new Size(200, 20),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Font = verdanaFont
            };
            sliceSelector.Items.AddRange(new string[] { "Surface", "Bottom", "Longitudinal", "Buoyancy", "Velocity", "Vertical Mixing" });
            sliceSelector.SelectedIndex = 0;
            sliceSelector.SelectedIndexChanged += SliceSelector_SelectedIndexChanged;

            // Buttons
            startButton = new Button
            {
                Location = new Point(10, 210),
                Size = new Size(200, 25),
                Text = "Start",
                FlatStyle = FlatStyle.Flat,
                Font = verdanaFont
            };
            startButton.Click += StartButton_Click;
            pauseButton = new Button
            {
                Location = new Point(10, 240),
                Size = new Size(200, 25),
                Text = "Pause",
                FlatStyle = FlatStyle.Flat,
                Enabled = false,
                Font = verdanaFont
            };
            pauseButton.Click += PauseButton_Click;
            resetButton = new Button
            {
                Location = new Point(10, 270),
                Size = new Size(200, 25),
                Text = "Reset",
                FlatStyle = FlatStyle.Flat,
                Font = verdanaFont
            };
            resetButton.Click += ResetButton_Click;

            // Visualization area
            visualizationBox = new PictureBox
            {
                Location = new Point(270, 10),
                Size = new Size(500, 400),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            visualizationBox.Paint += VisualizationBox_Paint;

            // Output console
            outputTextBox = new TextBox
            {
                Location = new Point(270, 420),
                Size = new Size(500, 130),
                Multiline = true,
                ReadOnly = true,
                ScrollBars = ScrollBars.Vertical,
                Font = verdanaFont
            };

            // Add controls to panel and form
            controlPanel.Controls.AddRange(new Control[] { mixingLabel, mixingCoefficientTextBox, tidalLabel, tidalAmplitudeTextBox, depthLabel, estuaryDepthTextBox, sliceLabel, sliceSelector, startButton, pauseButton, resetButton });
            turbulenceWindow.Controls.AddRange(new Control[] { controlPanel, visualizationBox, outputTextBox });

            // Center form on screen
            turbulenceWindow.StartPosition = FormStartPosition.CenterScreen;
        }

        private void InitializeSimulation()
        {
            salinity = new double[gridPointsX, gridPointsZ];
            u = new double[gridPointsX, gridPointsZ];
            w = new double[gridPointsX, gridPointsZ];
            p_h = new double[gridPointsX, gridPointsZ]; // Hydrostatic pressure
            p_nh = new double[gridPointsX, gridPointsZ]; // Non-hydrostatic pressure
            k = new double[gridPointsX, gridPointsZ];
            epsilon = new double[gridPointsX, gridPointsZ];
            buoyancy = new double[gridPointsX, gridPointsZ];
            K_x = new double[gridPointsX, gridPointsZ]; // Horizontal mixing coefficient
            K_z = new double[gridPointsX, gridPointsZ]; // Vertical mixing coefficient
            time = 0.0;
            isRunning = false;
            mixingCoefficient = 0.01;
            tidalAmplitude = 2.0;
            estuaryDepth = 10.0;
            dx = estuaryLength / gridPointsX;
            dz = estuaryDepth / gridPointsZ;
            dt = ComputeTimeStep();

            // Initialize fields with realistic estuarine profile
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    double x = i * dx;
                    double z = j * dz;
                    // Exponential salinity gradient (fresher upstream, saltier downstream)
                    double baseSalinity = maxSalinity * Math.Exp(-x / (estuaryLength / 2)) * (1.0 - z / estuaryDepth);
                    double perturbation = 0.05 * maxSalinity * (random.NextDouble() - 0.5); // ±2.5% perturbation
                    salinity[i, j] = Math.Max(0.0, Math.Min(maxSalinity, baseSalinity + perturbation));
                    // Logarithmic velocity profile near bottom
                    double frictionVelocity = 0.05; // Friction velocity (m/s)
                    u[i, j] = (z > z0) ? (frictionVelocity / vonKarman) * Math.Log(z / z0) : 0.0;
                    w[i, j] = 0.0;
                    // Initialize hydrostatic pressure by integrating buoyancy from surface downward
                    double rho = rho_0 + beta * (salinity[i, j] - refSalinity);
                    p_h[i, j] = j == gridPointsZ - 1 ? 0.0 : p_h[i, j + 1] + rho * g * dz;
                    p_nh[i, j] = 0.0; // Non-hydrostatic pressure starts at zero
                    k[i, j] = minK * (1 + 0.5 * Math.Exp(-z / estuaryDepth)); // Increased initial TKE
                    epsilon[i, j] = minEpsilon * (1 + 0.5 * Math.Exp(-z / estuaryDepth));
                    buoyancy[i, j] = -g * beta * (salinity[i, j] - refSalinity) / rho_0;
                    // Initialize anisotropic mixing coefficients
                    double eddyViscosity = C_mu * k[i, j] * k[i, j] / Math.Max(epsilon[i, j], minEpsilon);
                    K_x[i, j] = Math.Max(minMixing, Math.Min(maxMixing, eddyViscosity + mixingCoefficient));
                    K_z[i, j] = Math.Max(minMixing, Math.Min(maxMixing, eddyViscosity + mixingCoefficient));
                }
            }
            simulationTimer = new Timer
            {
                Interval = 100 // Update every 100ms
            };
            simulationTimer.Tick += (s, e) => UpdateSimulation();
        }

        private double ComputeTimeStep()
        {
            double maxU = 0.0, maxW = 0.0, maxK_x = 0.0, maxK_z = 0.0;
            for (int i = 0; i < gridPointsX; i++)
                for (int j = 0; j < gridPointsZ; j++)
                {
                    maxU = Math.Max(maxU, Math.Abs(u[i, j]));
                    maxW = Math.Max(maxW, Math.Abs(w[i, j]));
                    maxK_x = Math.Max(maxK_x, K_x[i, j]);
                    maxK_z = Math.Max(maxK_z, K_z[i, j]);
                }
            double cflAdv = cflFactor * Math.Min(dx / (maxU + 1e-6), dz / (maxW + 1e-6));
            double cflDiff = cflFactor * Math.Min(dx * dx / (2 * maxK_x + 1e-6), dz * dz / (2 * maxK_z + 1e-6));
            return Math.Min(0.02, Math.Min(cflAdv, cflDiff)); // Increased max dt for efficiency
        }

        private void UpdateSimulation()
        {
            try
            {
                // Parse and validate inputs
                mixingCoefficient = double.Parse(mixingCoefficientTextBox.Text);
                tidalAmplitude = double.Parse(tidalAmplitudeTextBox.Text);
                estuaryDepth = double.Parse(estuaryDepthTextBox.Text);
                if (mixingCoefficient < 1e-4 || mixingCoefficient > 0.1)
                    throw new Exception("Background mixing coefficient must be between 1e-4 and 0.1 m²/s.");
                if (tidalAmplitude < 0 || tidalAmplitude > 5.0)
                    throw new Exception("Tidal amplitude must be between 0 and 5 m.");
                if (estuaryDepth < 1.0 || estuaryDepth > 50.0)
                    throw new Exception("Estuary depth must be between 1 and 50 m.");

                dz = estuaryDepth / gridPointsZ;
                dt = ComputeTimeStep();
                double tidalVelocity = tidalAmplitude * Math.Sin(2 * Math.PI * time / 43200.0); // 12-hour tidal period

                // Temporary arrays
                double[,] uStar = new double[gridPointsX, gridPointsZ];
                double[,] wStar = new double[gridPointsX, gridPointsZ];
                double[,] newSalinity = new double[gridPointsX, gridPointsZ];
                double[,] newK = new double[gridPointsX, gridPointsZ];
                double[,] newEpsilon = new double[gridPointsX, gridPointsZ];
                double[,] newBuoyancy = new double[gridPointsX, gridPointsZ];
                double[,] newP_h = new double[gridPointsX, gridPointsZ];
                double[,] newK_x = new double[gridPointsX, gridPointsZ];
                double[,] newK_z = new double[gridPointsX, gridPointsZ];

                // Step 1: Compute intermediate velocity (uStar, wStar) without pressure and k-ε model
                for (int i = 1; i < gridPointsX - 1; i++)
                {
                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        k[i, j] = Math.Max(minK, k[i, j]);
                        epsilon[i, j] = Math.Max(minEpsilon, epsilon[i, j]);
                        // Compute eddy viscosity and anisotropic mixing coefficients
                        double eddyViscosity = C_mu * k[i, j] * k[i, j] / Math.Max(epsilon[i, j], minEpsilon);
                        // Gradient Richardson number for stratification effects
                        double ds_dz = j < gridPointsZ - 1 && j > 0 ? (salinity[i, j + 1] - salinity[i, j - 1]) / (2 * dz) : 0.0;
                        double N2 = -g * beta * ds_dz / rho_0; // Buoyancy frequency squared
                        // Compute velocity gradients for both momentum and k-ε model
                        double du_dx = i < gridPointsX - 1 && i > 0 ? (u[i + 1, j] - u[i - 1, j]) / (2 * dx) : 0.0;
                        double du_dz = j < gridPointsZ - 1 && j > 0 ? (u[i, j + 1] - u[i, j - 1]) / (2 * dz) : 0.0;
                        double dw_dx = i < gridPointsX - 1 && i > 0 ? (w[i + 1, j] - w[i - 1, j]) / (2 * dx) : 0.0;
                        double dw_dz = j < gridPointsZ - 1 && j > 0 ? (w[i, j + 1] - w[i, j - 1]) / (2 * dz) : 0.0;
                        double shear = du_dz * du_dz + dw_dz * dw_dz; // Simplified shear for Ri
                        double Ri = N2 / (shear + 1e-6); // Richardson number
                        double mixingReduction = Math.Max(0.1, Math.Min(1.0, 1.0 / (1.0 + 5.0 * Ri))); // Reduce vertical mixing
                        // Horizontal shear for K_x
                        double horizontalShear = du_dx * du_dx;
                        double shearEnhancement = Math.Min(2.0, 1.0 + 0.5 * horizontalShear / (1e-6 + horizontalShear));
                        newK_x[i, j] = Math.Max(minMixing, Math.Min(maxMixing, (eddyViscosity + mixingCoefficient) * shearEnhancement));
                        newK_z[i, j] = Math.Max(minMixing, Math.Min(maxMixing, (eddyViscosity + mixingCoefficient) * mixingReduction));

                        // Buoyancy
                        double salinityPerturbation = salinity[i, j] - refSalinity;
                        double buoyancyValue = Math.Max(-5e-2, Math.Min(5e-2, -g * beta * salinityPerturbation / rho_0));
                        newBuoyancy[i, j] = buoyancyValue;

                        // Update hydrostatic pressure
                        double rho = rho_0 + beta * (salinity[i, j] - refSalinity);
                        newP_h[i, j] = j == gridPointsZ - 1 ? 0.0 : newP_h[i, j + 1] + rho * g * dz;

                        // Upwind advection terms with boundary handling
                        double advection_u_x, advection_u_z, advection_w_x, advection_w_z;
                        if (i >= 2 && i < gridPointsX - 2)
                        {
                            advection_u_x = u[i, j] > 0 ? (3 * u[i, j] - 4 * u[i - 1, j] + u[i - 2, j]) / (2 * dx) :
                                                          (u[i + 2, j] - 4 * u[i + 1, j] + 3 * u[i, j]) / (2 * dx);
                            advection_w_x = u[i, j] > 0 ? (3 * w[i, j] - 4 * w[i - 1, j] + w[i - 2, j]) / (2 * dx) :
                                                          (w[i + 2, j] - 4 * w[i + 1, j] + 3 * w[i, j]) / (2 * dx);
                        }
                        else
                        {
                            advection_u_x = u[i, j] > 0 ? (u[i, j] - u[i - 1, j]) / dx : (u[i + 1, j] - u[i, j]) / dx;
                            advection_w_x = u[i, j] > 0 ? (w[i, j] - w[i - 1, j]) / dx : (w[i + 1, j] - w[i, j]) / dx;
                        }
                        if (j >= 2 && j < gridPointsZ - 2)
                        {
                            advection_u_z = w[i, j] > 0 ? (3 * u[i, j] - 4 * u[i, j - 1] + u[i, j - 2]) / (2 * dz) :
                                                          (u[i, j + 2] - 4 * u[i, j + 1] + 3 * u[i, j]) / (2 * dz);
                            advection_w_z = w[i, j] > 0 ? (3 * w[i, j] - 4 * w[i, j - 1] + w[i, j - 2]) / (2 * dz) :
                                                          (w[i, j + 2] - 4 * w[i, j + 1] + 3 * w[i, j]) / (2 * dz);
                        }
                        else
                        {
                            advection_u_z = w[i, j] > 0 ? (u[i, j] - u[i, j - 1]) / dz : (u[i, j + 1] - u[i, j]) / dz;
                            advection_w_z = w[i, j] > 0 ? (w[i, j] - w[i, j - 1]) / dz : (w[i, j + 1] - w[i, j]) / dz;
                        }

                        // Anisotropic diffusion terms
                        double diffusion_u = newK_x[i, j] * (u[i + 1, j] - 2 * u[i, j] + u[i - 1, j]) / (dx * dx) +
                                            newK_z[i, j] * (u[i, j + 1] - 2 * u[i, j] + u[i, j - 1]) / (dz * dz);
                        double diffusion_w = newK_x[i, j] * (w[i + 1, j] - 2 * w[i, j] + w[i - 1, j]) / (dx * dx) +
                                            newK_z[i, j] * (w[i, j + 1] - 2 * w[i, j] + w[i, j - 1]) / (dz * dz);

                        uStar[i, j] = Math.Max(-2.0, Math.Min(2.0, u[i, j] + dt * (-advection_u_x - advection_u_z + diffusion_u)));
                        wStar[i, j] = Math.Max(-0.2, Math.Min(0.2, w[i, j] + dt * (-advection_w_x - advection_w_z + diffusion_w + buoyancyValue)));
                    }
                }

                // Step 2: Solve non-hydrostatic pressure Poisson equation
                double[,] newP_nh = new double[gridPointsX, gridPointsZ];
                for (int iter = 0; iter < 100; iter++)
                {
                    for (int i = 1; i < gridPointsX - 1; i++)
                    {
                        for (int j = 1; j < gridPointsZ - 1; j++)
                        {
                            double divU = (uStar[i + 1, j] - uStar[i - 1, j]) / (2 * dx) +
                                          (wStar[i, j + 1] - wStar[i, j - 1]) / (2 * dz);
                            double laplacianP = (p_nh[i + 1, j] + p_nh[i - 1, j]) / (dx * dx) +
                                                (p_nh[i, j + 1] + p_nh[i, j - 1]) / (dz * dz) -
                                                2 * p_nh[i, j] * (1 / (dx * dx) + 1 / (dz * dz));
                            newP_nh[i, j] = p_nh[i, j] + pressureRelaxation * (laplacianP - rho_0 * divU / dt);
                            newP_nh[i, j] = Math.Max(-1e4, Math.Min(1e4, newP_nh[i, j]));
                        }
                    }
                    p_nh = newP_nh;
                }

                // Step 3: Correct velocity with total pressure (hydrostatic + non-hydrostatic)
                for (int i = 1; i < gridPointsX - 1; i++)
                {
                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        double dp_dx = (p_h[i + 1, j] + p_nh[i + 1, j] - p_h[i - 1, j] - p_nh[i - 1, j]) / (2 * dx);
                        double dp_dz = (p_h[i, j + 1] + p_nh[i, j + 1] - p_h[i, j - 1] - p_nh[i, j - 1]) / (2 * dz);
                        u[i, j] = Math.Max(-2.0, Math.Min(2.0, uStar[i, j] - dt / rho_0 * dp_dx));
                        w[i, j] = Math.Max(-0.2, Math.Min(0.2, wStar[i, j] - dt / rho_0 * dp_dz));
                        u[i, j] += dt * (tidalVelocity - u[i, j]);
                    }
                }

                // Step 4: Update k-ε model and salinity
                for (int i = 1; i < gridPointsX - 1; i++)
                {
                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        // Shear production (reusing velocity gradients)
                        double du_dx = (u[i + 1, j] - u[i - 1, j]) / (2 * dx);
                        double du_dz = (u[i, j + 1] - u[i, j - 1]) / (2 * dz);
                        double dw_dx = (w[i + 1, j] - w[i - 1, j]) / (2 * dx);
                        double dw_dz = (w[i, j + 1] - w[i, j - 1]) / (2 * dz);
                        double production = Math.Min(0.1, (newK_x[i, j] + newK_z[i, j]) / 2.0 *
                            (2 * du_dx * du_dx + 2 * dw_dz * dw_dz + (du_dz + dw_dx) * (du_dz + dw_dx))); // Full shear tensor

                        // Buoyancy
                        double salinityPerturbation = salinity[i, j] - refSalinity;
                        double buoyancyValue = Math.Max(-5e-2, Math.Min(5e-2, -g * beta * salinityPerturbation / rho_0));
                        newBuoyancy[i, j] = buoyancyValue;

                        // k-equation with anisotropic diffusion
                        double diffusion_k = (newK_x[i, j] / sigma_k) * (k[i + 1, j] - 2 * k[i, j] + k[i - 1, j]) / (dx * dx) +
                                            (newK_z[i, j] / sigma_k) * (k[i, j + 1] - 2 * k[i, j] + k[i, j - 1]) / (dz * dz);
                        newK[i, j] = Math.Max(minK, k[i, j] + dt * (production - epsilon[i, j] + diffusion_k - buoyancyValue));

                        // ε-equation with anisotropic diffusion
                        double diffusion_epsilon = (newK_x[i, j] / sigma_epsilon) * (epsilon[i + 1, j] - 2 * epsilon[i, j] + epsilon[i - 1, j]) / (dx * dx) +
                                                 (newK_z[i, j] / sigma_epsilon) * (epsilon[i, j + 1] - 2 * epsilon[i, j] + epsilon[i, j - 1]) / (dz * dz);
                        double epsilonTerm = C1_epsilon * production * epsilon[i, j] / Math.Max(k[i, j], minK) -
                                            C2_epsilon * epsilon[i, j] * epsilon[i, j] / Math.Max(k[i, j], minK);
                        newEpsilon[i, j] = Math.Max(minEpsilon, epsilon[i, j] + dt * (epsilonTerm + diffusion_epsilon));

                        // Salinity transport with anisotropic diffusion
                        double advection_s_x, advection_s_z;
                        if (i >= 2 && i < gridPointsX - 2)
                        {
                            advection_s_x = u[i, j] > 0 ? (3 * salinity[i, j] - 4 * salinity[i - 1, j] + salinity[i - 2, j]) / (2 * dx) :
                                                          (salinity[i + 2, j] - 4 * salinity[i + 1, j] + 3 * salinity[i, j]) / (2 * dx);
                        }
                        else
                        {
                            advection_s_x = u[i, j] > 0 ? (salinity[i, j] - salinity[i - 1, j]) / dx :
                                                          (salinity[i + 1, j] - salinity[i, j]) / dx;
                        }
                        if (j >= 2 && j < gridPointsZ - 2)
                        {
                            advection_s_z = w[i, j] > 0 ? (3 * salinity[i, j] - 4 * salinity[i, j - 1] + salinity[i, j - 2]) / (2 * dz) :
                                                          (salinity[i, j + 2] - 4 * salinity[i, j + 1] + 3 * salinity[i, j]) / (2 * dz);
                        }
                        else
                        {
                            advection_s_z = w[i, j] > 0 ? (salinity[i, j] - salinity[i, j - 1]) / dz :
                                                          (salinity[i, j + 1] - salinity[i, j]) / dz;
                        }
                        double diffusion_s = newK_x[i, j] * (salinity[i + 1, j] - 2 * salinity[i, j] + salinity[i - 1, j]) / (dx * dx) +
                                            newK_z[i, j] * (salinity[i, j + 1] - 2 * salinity[i, j] + salinity[i, j - 1]) / (dz * dz);
                        newSalinity[i, j] = Math.Max(0.0, Math.Min(maxSalinity, salinity[i, j] + dt * (diffusion_s - advection_s_x - advection_s_z)));
                    }
                }

                // Boundary conditions
                for (int i = 0; i < gridPointsX; i++)
                {
                    double frictionVelocity = 0.05; // Friction velocity
                    newK[i, 0] = newK[i, 1];
                    newK[i, gridPointsZ - 1] = newK[i, gridPointsZ - 2];
                    newEpsilon[i, 0] = Math.Pow(C_mu, 0.75) * Math.Pow(newK[i, 1], 1.5) / (vonKarman * dz); // Improved wall function
                    newEpsilon[i, gridPointsZ - 1] = newEpsilon[i, gridPointsZ - 2];
                    newSalinity[i, 0] = newSalinity[i, 1];
                    newSalinity[i, gridPointsZ - 1] = maxSalinity;
                    newBuoyancy[i, 0] = newBuoyancy[i, 1];
                    newBuoyancy[i, gridPointsZ - 1] = -g * beta * (newSalinity[i, gridPointsZ - 1] - refSalinity) / rho_0;
                    newP_h[i, 0] = newP_h[i, 1] + (rho_0 + beta * (newSalinity[i, 0] - refSalinity)) * g * dz;
                    newP_h[i, gridPointsZ - 1] = 0.0; // Surface pressure is zero
                    p_nh[i, 0] = p_nh[i, 1]; // Neumann condition for non-hydrostatic pressure
                    p_nh[i, gridPointsZ - 1] = 0.0; // Zero non-hydrostatic pressure at surface
                    newK_x[i, 0] = newK_x[i, 1] * 0.1; // Reduced mixing near bottom
                    newK_x[i, gridPointsZ - 1] = newK_x[i, gridPointsZ - 2];
                    newK_z[i, 0] = newK_z[i, 1] * 0.01; // Strongly reduced vertical mixing near bottom
                    newK_z[i, gridPointsZ - 1] = newK_z[i, gridPointsZ - 2];
                    u[i, 0] = (frictionVelocity / vonKarman) * Math.Log(dz / z0); // Logarithmic profile
                    u[i, gridPointsZ - 1] = tidalVelocity;
                    w[i, 0] = 0.0;
                    w[i, gridPointsZ - 1] = 0.0;
                }
                for (int j = 0; j < gridPointsZ; j++)
                {
                    newK[0, j] = newK[1, j];
                    newK[gridPointsX - 1, j] = newK[gridPointsX - 2, j];
                    newEpsilon[0, j] = newEpsilon[1, j];
                    newEpsilon[gridPointsX - 1, j] = newEpsilon[gridPointsX - 2, j];
                    newSalinity[0, j] = 0.0;
                    newSalinity[gridPointsX - 1, j] = newSalinity[gridPointsX - 2, j];
                    newBuoyancy[0, j] = -g * beta * (newSalinity[0, j] - refSalinity) / rho_0;
                    newBuoyancy[gridPointsX - 1, j] = newBuoyancy[gridPointsX - 2, j];
                    newP_h[0, j] = newP_h[1, j]; // Neumann condition for hydrostatic pressure
                    newP_h[gridPointsX - 1, j] = newP_h[gridPointsX - 2, j];
                    p_nh[0, j] = p_nh[1, j]; // Neumann condition for non-hydrostatic pressure
                    p_nh[gridPointsX - 1, j] = p_nh[gridPointsX - 2, j];
                    newK_x[0, j] = newK_x[1, j];
                    newK_x[gridPointsX - 1, j] = newK_x[gridPointsX - 2, j];
                    newK_z[0, j] = newK_z[1, j];
                    newK_z[gridPointsX - 1, j] = newK_z[gridPointsX - 2, j];
                    u[0, j] = 0.05;
                    u[gridPointsX - 1, j] = tidalVelocity;
                    w[0, j] = 0.0;
                    w[gridPointsX - 1, j] = 0.0;
                }

                k = newK;
                epsilon = newEpsilon;
                salinity = newSalinity;
                buoyancy = newBuoyancy;
                p_h = newP_h;
                K_x = newK_x;
                K_z = newK_z;

                time += dt;
                visualizationBox.Invalidate();
                UpdateOutputConsole();
            }
            catch (Exception ex)
            {
                simulationTimer.Stop();
                isRunning = false;
                MessageBox.Show($"Simulation Error: {ex.Message}", "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                UpdateButtonStates();
            }
        }

        private void UpdateOutputConsole()
        {
            double avgSalinity = 0.0, avgK = 0.0, avgEpsilon = 0.0, avgUSpeed = 0.0, avgK_x = 0.0, avgK_z = 0.0;
            for (int i = 0; i < gridPointsX; i++)
                for (int j = 0; j < gridPointsZ; j++)
                {
                    avgSalinity += salinity[i, j];
                    avgK += k[i, j];
                    avgEpsilon += epsilon[i, j];
                    avgUSpeed += Math.Sqrt(u[i, j] * u[i, j] + w[i, j] * w[i, j]);
                    avgK_x += K_x[i, j];
                    avgK_z += K_z[i, j];
                }
            avgSalinity /= (gridPointsX * gridPointsZ);
            avgK /= (gridPointsX * gridPointsZ);
            avgEpsilon /= (gridPointsX * gridPointsZ);
            avgUSpeed /= (gridPointsX * gridPointsZ);
            avgK_x /= (gridPointsX * gridPointsZ);
            avgK_z /= (gridPointsX * gridPointsZ);

            outputTextBox.AppendText($"Time: {time:F2}s | Avg Salinity: {avgSalinity:F2} PSU | Avg TKE: {avgK:E2} m²/s² | Avg Epsilon: {avgEpsilon:E2} m²/s³ | Avg Speed: {avgUSpeed:F4} m/s | Avg K_x: {avgK_x:F4} m²/s | Avg K_z: {avgK_z:F4} m²/s | dt: {dt:F6}s\r\n");
        }

        private void VisualizationBox_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            int width = visualizationBox.Width;
            int height = visualizationBox.Height;
            double minSalinity = 0.0;

            if (currentSlice == "Surface" || currentSlice == "Bottom")
            {
                int zIndex = currentSlice == "Surface" ? gridPointsZ - 1 : 0;
                for (int i = 0; i < gridPointsX; i++)
                {
                    for (int px = i * width / gridPointsX; px < (i + 1) * width / gridPointsX; px++)
                    {
                        for (int py = 0; py < height; py++)
                        {
                            float normalized = (float)((salinity[i, zIndex] - minSalinity) / (maxSalinity - minSalinity));
                            normalized = Math.Max(0.0f, Math.Min(1.0f, normalized));
                            Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                            g.FillRectangle(new SolidBrush(color), px, py, 1, 1);
                        }
                    }
                }
                // Enhanced color bar
                for (int y = 0; y < height; y++)
                {
                    float normalized = 1.0f - (float)y / height;
                    Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                    g.FillRectangle(new SolidBrush(color), width - 20, y, 20, 1);
                }
                g.DrawString("35 PSU", verdanaFont, Brushes.Black, width - 40, 0);
                g.DrawString("0 PSU", verdanaFont, Brushes.Black, width - 40, height - 20);
            }
            else if (currentSlice == "Buoyancy")
            {
                int zIndex = gridPointsZ - 1;
                for (int i = 0; i < gridPointsX; i++)
                {
                    for (int px = i * width / gridPointsX; px < (i + 1) * width / gridPointsX; px++)
                    {
                        for (int py = 0; py < height; py++)
                        {
                            float normalized = (float)((buoyancy[i, zIndex] - minBuoyancy) / (maxBuoyancy - minBuoyancy));
                            normalized = Math.Max(0.0f, Math.Min(1.0f, normalized));
                            Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                            g.FillRectangle(new SolidBrush(color), px, py, 1, 1);
                        }
                    }
                }
                for (int y = 0; y < height; y++)
                {
                    float normalized = 1.0f - (float)y / height;
                    Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                    g.FillRectangle(new SolidBrush(color), width - 20, y, 20, 1);
                }
                g.DrawString("0 m/s²", verdanaFont, Brushes.Black, width - 40, 0);
                g.DrawString("-0.3 m/s²", verdanaFont, Brushes.Black, width - 40, height - 20);
            }
            else if (currentSlice == "Velocity")
            {
                // Visualize velocity magnitude with vectors
                double maxSpeed = 0.0;
                for (int i = 0; i < gridPointsX; i++)
                    for (int j = 0; j < gridPointsZ; j++)
                        maxSpeed = Math.Max(maxSpeed, Math.Sqrt(u[i, j] * u[i, j] + w[i, j] * w[i, j]));
                for (int i = 0; i < gridPointsX; i += 5)
                {
                    for (int j = 0; j < gridPointsZ; j += 5)
                    {
                        int px = i * width / gridPointsX;
                        int py = height - (j * height / gridPointsZ);
                        float speed = (float)Math.Sqrt(u[i, j] * u[i, j] + w[i, j] * w[i, j]);
                        float normalized = (float)(speed / (maxSpeed + 1e-6));
                        Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                        g.FillRectangle(new SolidBrush(color), px, py, width / gridPointsX, height / gridPointsZ);
                        // Draw velocity vector
                        float scale = 20.0f;
                        int endX = px + (int)(u[i, j] * scale);
                        int endY = py - (int)(w[i, j] * scale);
                        g.DrawLine(Pens.Black, px, py, endX, endY);
                    }
                }
                for (int y = 0; y < height; y++)
                {
                    float normalized = 1.0f - (float)y / height;
                    Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                    g.FillRectangle(new SolidBrush(color), width - 20, y, 20, 1);
                }
                g.DrawString($"{maxSpeed:F2} m/s", verdanaFont, Brushes.Black, width - 40, 0);
                g.DrawString("0 m/s", verdanaFont, Brushes.Black, width - 40, height - 20);
            }
            else if (currentSlice == "VerticalMixing")
            {
                // Visualize K_z (vertical mixing coefficient)
                double minK_z = minMixing;
                double maxK_z = maxMixing;
                for (int i = 0; i < gridPointsX; i++)
                {
                    for (int j = 0; j < gridPointsZ; j++)
                    {
                        int px = i * width / gridPointsX;
                        int py = height - (j * height / gridPointsZ);
                        float normalized = (float)((K_z[i, j] - minK_z) / (maxK_z - minK_z));
                        normalized = Math.Max(0.0f, Math.Min(1.0f, normalized));
                        Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                        g.FillRectangle(new SolidBrush(color), px, py, width / gridPointsX, height / gridPointsZ);
                    }
                }
                for (int y = 0; y < height; y++)
                {
                    float normalized = 1.0f - (float)y / height;
                    Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                    g.FillRectangle(new SolidBrush(color), width - 20, y, 20, 1);
                }
                g.DrawString($"{maxK_z:F4} m²/s", verdanaFont, Brushes.Black, width - 40, 0);
                g.DrawString($"{minK_z:F4} m²/s", verdanaFont, Brushes.Black, width - 40, height - 20);
            }
            else // Longitudinal
            {
                for (int i = 0; i < gridPointsX; i++)
                {
                    for (int j = 0; j < gridPointsZ; j++)
                    {
                        int px = i * width / gridPointsX;
                        int py = height - (j * height / gridPointsZ);
                        float normalized = (float)((salinity[i, j] - minSalinity) / (maxSalinity - minSalinity));
                        normalized = Math.Max(0.0f, Math.Min(1.0f, normalized));
                        Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                        g.FillRectangle(new SolidBrush(color), px, py, width / gridPointsX, height / gridPointsZ);
                    }
                }
                for (int y = 0; y < height; y++)
                {
                    float normalized = 1.0f - (float)y / height;
                    Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                    g.FillRectangle(new SolidBrush(color), width - 20, y, 20, 1);
                }
                g.DrawString("35 PSU", verdanaFont, Brushes.Black, width - 40, 0);
                g.DrawString("0 PSU", verdanaFont, Brushes.Black, width - 40, height - 20);
            }
        }

        private void SliceSelector_SelectedIndexChanged(object sender, EventArgs e)
        {
            currentSlice = sliceSelector.SelectedItem.ToString();
            visualizationBox.Invalidate();
        }

        private void StartButton_Click(object sender, EventArgs e)
        {
            try
            {
                mixingCoefficient = double.Parse(mixingCoefficientTextBox.Text);
                tidalAmplitude = double.Parse(tidalAmplitudeTextBox.Text);
                estuaryDepth = double.Parse(estuaryDepthTextBox.Text);
                if (mixingCoefficient < 1e-4 || mixingCoefficient > 0.1)
                    throw new Exception("Background mixing coefficient must be between 1e-4 and 0.1 m²/s.");
                if (tidalAmplitude < 0 || tidalAmplitude > 5.0)
                    throw new Exception("Tidal amplitude must be between 0 and 5 m.");
                if (estuaryDepth < 1.0 || estuaryDepth > 50.0)
                    throw new Exception("Estuary depth must be between 1 and 50 m.");

                dz = estuaryDepth / gridPointsZ;
                dt = ComputeTimeStep();
                simulationTimer.Start();
                isRunning = true;
                outputTextBox.AppendText("Turbulence and Mixing simulation started.\r\n");
                UpdateButtonStates();
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
            outputTextBox.AppendText("Simulation paused.\r\n");
            UpdateButtonStates();
        }

        private void ResetButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isRunning = false;
            InitializeSimulation();
            outputTextBox.Clear();
            outputTextBox.AppendText("Simulation reset.\r\n");
            visualizationBox.Invalidate();
            UpdateButtonStates();
        }

        private void UpdateButtonStates()
        {
            startButton.Enabled = !isRunning;
            pauseButton.Enabled = isRunning;
            resetButton.Enabled = true;
        }

        public void ShowTurbulenceWindow()
        {
            turbulenceWindow.Show();
        }
    }
}