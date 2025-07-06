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
        private CheckBox useThirdOrderAdvectionCheckBox; // New checkbox
        private Timer simulationTimer;
        private double[,] salinity; // 2D salinity field (ξ, η)
        private double[,] u; // Horizontal velocity (m/s)
        private double[,] w; // Vertical velocity (m/s)
        private double[,] p_h; // Hydrostatic pressure (Pa)
        private double[,] p_nh; // Non-hydrostatic pressure (Pa)
        private double[,] k; // Turbulent kinetic energy
        private double[,] epsilon; // Dissipation rate
        private double[,] buoyancy; // Buoyancy field (m/s²)
        private double[,] K_x; // Horizontal turbulent mixing coefficient (m²/s)
        private double[,] K_z; // Vertical turbulent mixing coefficient (m²/s)
        private double[,] x; // Physical x-coordinate (m)
        private double[,] z; // Physical z-coordinate (m)
        private double[,] J; // Jacobian of coordinate transformation
        private double[,] x_xi; // ∂x/∂ξ
        private double[,] x_eta; // ∂x/∂η
        private double[,] z_xi; // ∂z/∂ξ
        private double[,] z_eta; // ∂z/∂η
        private double mixingCoefficient; // Background mixing coefficient
        private double tidalAmplitude; // Tidal amplitude (m)
        private double estuaryDepth; // Estuary depth (m)
        private double estuaryLength = 10000.0; // Fixed estuary length (m)
        private int gridPointsX = 100; // Logical grid points in ξ
        private int gridPointsZ = 50; // Logical grid points in η
        private double dxi; // Logical grid spacing in ξ
        private double deta; // Logical grid spacing in η
        private double time; // Simulation time
        private double dt; // Adaptive time step
        private bool isRunning;
        private string currentSlice = "Surface"; // Default visualization
        private bool useThirdOrderAdvection = true; // Default to third-order scheme
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
        private const double vonKarman = 0.41; // von Karman constant
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
            sliceSelector.Items.AddRange(new string[] { "Surface", "Bottom", "Longitudinal", "Buoyancy", "Velocity", "VerticalMixing" });
            sliceSelector.SelectedIndex = 0;
            sliceSelector.Items.Add("Advection");
            sliceSelector.SelectedIndexChanged += SliceSelector_SelectedIndexChanged;

            // Checkbox for advection scheme
            useThirdOrderAdvectionCheckBox = new CheckBox
            {
                Location = new Point(10, 210),
                Size = new Size(200, 20),
                Text = "Use Third-Order Advection",
                Checked = true,
                Font = verdanaFont
            };
            useThirdOrderAdvectionCheckBox.CheckedChanged += UseThirdOrderAdvectionCheckBox_CheckedChanged;

            // Buttons
            startButton = new Button
            {
                Location = new Point(10, 240),
                Size = new Size(200, 25),
                Text = "Start",
                FlatStyle = FlatStyle.Flat,
                Font = verdanaFont
            };
            startButton.Click += StartButton_Click;
            pauseButton = new Button
            {
                Location = new Point(10, 270),
                Size = new Size(200, 25),
                Text = "Pause",
                FlatStyle = FlatStyle.Flat,
                Enabled = false,
                Font = verdanaFont
            };
            pauseButton.Click += PauseButton_Click;
            resetButton = new Button
            {
                Location = new Point(10, 300),
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
            controlPanel.Controls.AddRange(new Control[] { mixingLabel, mixingCoefficientTextBox, tidalLabel, tidalAmplitudeTextBox, depthLabel, estuaryDepthTextBox, sliceLabel, sliceSelector, useThirdOrderAdvectionCheckBox, startButton, pauseButton, resetButton });
            turbulenceWindow.Controls.AddRange(new Control[] { controlPanel, visualizationBox, outputTextBox });

            // Center form on screen
            turbulenceWindow.StartPosition = FormStartPosition.CenterScreen;
        }

        private void InitializeSimulation()
        {
            salinity = new double[gridPointsX, gridPointsZ];
            u = new double[gridPointsX, gridPointsZ];
            w = new double[gridPointsX, gridPointsZ];
            p_h = new double[gridPointsX, gridPointsZ];
            p_nh = new double[gridPointsX, gridPointsZ];
            k = new double[gridPointsX, gridPointsZ];
            epsilon = new double[gridPointsX, gridPointsZ];
            buoyancy = new double[gridPointsX, gridPointsZ];
            K_x = new double[gridPointsX, gridPointsZ];
            K_z = new double[gridPointsX, gridPointsZ];
            x = new double[gridPointsX, gridPointsZ];
            z = new double[gridPointsX, gridPointsZ];
            J = new double[gridPointsX, gridPointsZ];
            x_xi = new double[gridPointsX, gridPointsZ];
            x_eta = new double[gridPointsX, gridPointsZ];
            z_xi = new double[gridPointsX, gridPointsZ];
            z_eta = new double[gridPointsX, gridPointsZ];
            time = 0.0;
            isRunning = false;
            mixingCoefficient = 0.01;
            tidalAmplitude = 2.0;
            estuaryDepth = 10.0;
            dxi = 1.0 / (gridPointsX - 1); // Logical grid spacing in ξ
            deta = 1.0 / (gridPointsZ - 1); // Logical grid spacing in η
            dt = ComputeTimeStep();

            // Define curvilinear grid: x(ξ, η), z(ξ, η)
            double widthAtMouth = 2000.0; // Width at x = estuaryLength (m)
            double widthAtHead = 500.0; // Width at x = 0 (m)
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    double xi = i * dxi;
                    double eta = j * deta;
                    // Exponential width variation
                    double width = widthAtHead + (widthAtMouth - widthAtHead) * Math.Exp(-2.0 * (1.0 - xi));
                    x[i, j] = xi * estuaryLength + width * eta * (1.0 - eta); // Lateral distortion
                    z[i, j] = eta * estuaryDepth; // Depth follows η
                }
            }

            // Compute metric terms
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    double xi_x = i > 0 && i < gridPointsX - 1 ? (x[i + 1, j] - x[i - 1, j]) / (2 * dxi) :
                        i == 0 ? (x[1, j] - x[0, j]) / dxi : (x[gridPointsX - 1, j] - x[gridPointsX - 2, j]) / dxi;
                    double xi_z = i > 0 && i < gridPointsX - 1 ? (z[i + 1, j] - z[i - 1, j]) / (2 * dxi) :
                        i == 0 ? (z[1, j] - z[0, j]) / dxi : (z[gridPointsX - 1, j] - z[gridPointsX - 2, j]) / dxi;
                    double eta_x = j > 0 && j < gridPointsZ - 1 ? (x[i, j + 1] - x[i, j - 1]) / (2 * deta) :
                        j == 0 ? (x[i, 1] - x[i, 0]) / deta : (x[i, gridPointsZ - 1] - x[i, gridPointsZ - 2]) / deta;
                    double eta_z = j > 0 && j < gridPointsZ - 1 ? (z[i, j + 1] - z[i, j - 1]) / (2 * deta) :
                        j == 0 ? (z[i, 1] - z[i, 0]) / deta : (z[i, gridPointsZ - 1] - z[i, gridPointsZ - 2]) / deta;
                    x_xi[i, j] = xi_x;
                    x_eta[i, j] = eta_x;
                    z_xi[i, j] = xi_z;
                    z_eta[i, j] = eta_z;
                    J[i, j] = xi_x * eta_z - xi_z * eta_x; // Jacobian
                    if (Math.Abs(J[i, j]) < 1e-6) J[i, j] = 1e-6; // Prevent division by zero
                }
            }

            // Initialize fields
            for (int i = 0; i < gridPointsX; i++)
            {
                for (int j = 0; j < gridPointsZ; j++)
                {
                    double phys_x = x[i, j];
                    double phys_z = z[i, j];
                    // Exponential salinity gradient
                    double baseSalinity = maxSalinity * Math.Exp(-phys_x / (estuaryLength / 2)) * (1.0 - phys_z / estuaryDepth);
                    double perturbation = 0.05 * maxSalinity * (random.NextDouble() - 0.5);
                    salinity[i, j] = Math.Max(0.0, Math.Min(maxSalinity, baseSalinity + perturbation));
                    double frictionVelocity = 0.05;
                    u[i, j] = (phys_z > z0) ? (frictionVelocity / vonKarman) * Math.Log(phys_z / z0) : 0.0;
                    w[i, j] = 0.0;
                    double rho = rho_0 + beta * (salinity[i, j] - refSalinity);
                    p_h[i, j] = j == gridPointsZ - 1 ? 0.0 : p_h[i, j + 1] + rho * g * z_eta[i, j] * deta;
                    p_nh[i, j] = 0.0;
                    k[i, j] = minK * (1 + 0.5 * Math.Exp(-phys_z / estuaryDepth));
                    epsilon[i, j] = minEpsilon * (1 + 0.5 * Math.Exp(-phys_z / estuaryDepth));
                    buoyancy[i, j] = -g * beta * (salinity[i, j] - refSalinity) / rho_0;
                    double eddyViscosity = C_mu * k[i, j] * k[i, j] / Math.Max(epsilon[i, j], minEpsilon);
                    K_x[i, j] = Math.Max(minMixing, Math.Min(maxMixing, eddyViscosity + mixingCoefficient));
                    K_z[i, j] = Math.Max(minMixing, Math.Min(maxMixing, eddyViscosity + mixingCoefficient));
                }
            }
            simulationTimer = new Timer
            {
                Interval = 100
            };
            simulationTimer.Tick += (s, e) => UpdateSimulation();
        }

        private double ComputeTimeStep()
        {
            double maxU = 0.0, maxW = 0.0, maxK_x = 0.0, maxK_z = 0.0;
            double min_dx = double.MaxValue, min_dz = double.MaxValue;
            for (int i = 0; i < gridPointsX; i++)
                for (int j = 0; j < gridPointsZ; j++)
                {
                    maxU = Math.Max(maxU, Math.Abs(u[i, j]));
                    maxW = Math.Max(maxW, Math.Abs(w[i, j]));
                    maxK_x = Math.Max(maxK_x, K_x[i, j]);
                    maxK_z = Math.Max(maxK_z, K_z[i, j]);
                    if (i > 0) min_dx = Math.Min(min_dx, Math.Abs(x[i, j] - x[i - 1, j]));
                    if (j > 0) min_dz = Math.Min(min_dz, Math.Abs(z[i, j] - z[i, j - 1]));
                }
            double cflFactor = useThirdOrderAdvection ? 0.5 : 0.7; // Smaller CFL for third-order to ensure stability
            double cflAdv = cflFactor * Math.Min(min_dx / (maxU + 1e-6), min_dz / (maxW + 1e-6));
            double cflDiff = cflFactor * Math.Min(min_dx * min_dx / (2 * maxK_x + 1e-6), min_dz * min_dz / (2 * maxK_z + 1e-6));
            return Math.Min(0.02, Math.Min(cflAdv, cflDiff));
        }

        private void UpdateSimulation()
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

                deta = estuaryDepth / (gridPointsZ - 1);
                dt = ComputeTimeStep();
                double tidalVelocity = tidalAmplitude * Math.Sin(2 * Math.PI * time / 43200.0);

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

                // Step 1: Compute intermediate velocity
                for (int i = 1; i < gridPointsX - 1; i++)
                {
                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        k[i, j] = Math.Max(minK, k[i, j]);
                        epsilon[i, j] = Math.Max(minEpsilon, epsilon[i, j]);
                        double eddyViscosity = C_mu * k[i, j] * k[i, j] / Math.Max(epsilon[i, j], minEpsilon);
                        // Gradients in logical coordinates
                        double ds_dxi = (salinity[i + 1, j] - salinity[i - 1, j]) / (2 * dxi);
                        double ds_deta = (salinity[i, j + 1] - salinity[i, j - 1]) / (2 * deta);
                        double ds_dx = (ds_dxi * z_eta[i, j] - ds_deta * z_xi[i, j]) / J[i, j];
                        double ds_dz = (-ds_dxi * x_eta[i, j] + ds_deta * x_xi[i, j]) / J[i, j];
                        double N2 = -g * beta * ds_dz / rho_0;
                        double du_dxi = (u[i + 1, j] - u[i - 1, j]) / (2 * dxi);
                        double du_deta = (u[i, j + 1] - u[i, j - 1]) / (2 * deta);
                        double dw_dxi = (w[i + 1, j] - w[i - 1, j]) / (2 * dxi);
                        double dw_deta = (w[i, j + 1] - w[i, j - 1]) / (2 * deta);
                        double du_dx = (du_dxi * z_eta[i, j] - du_deta * z_xi[i, j]) / J[i, j];
                        double du_dz = (-du_dxi * x_eta[i, j] + du_deta * x_xi[i, j]) / J[i, j];
                        double dw_dx = (dw_dxi * z_eta[i, j] - dw_deta * z_xi[i, j]) / J[i, j];
                        double dw_dz = (-dw_dxi * x_eta[i, j] + dw_deta * x_xi[i, j]) / J[i, j];
                        double shear = du_dz * du_dz + dw_dz * dw_dz;
                        double Ri = N2 / (shear + 1e-6);
                        double mixingReduction = Math.Max(0.1, Math.Min(1.0, 1.0 / (1.0 + 5.0 * Ri)));
                        double horizontalShear = du_dx * du_dx;
                        double shearEnhancement = Math.Min(2.0, 1.0 + 0.5 * horizontalShear / (1e-6 + horizontalShear));
                        newK_x[i, j] = Math.Max(minMixing, Math.Min(maxMixing, (eddyViscosity + mixingCoefficient) * shearEnhancement));
                        newK_z[i, j] = Math.Max(minMixing, Math.Min(maxMixing, (eddyViscosity + mixingCoefficient) * mixingReduction));

                        double salinityPerturbation = salinity[i, j] - refSalinity;
                        double buoyancyValue = Math.Max(-5e-2, Math.Min(5e-2, -g * beta * salinityPerturbation / rho_0));
                        newBuoyancy[i, j] = buoyancyValue;

                        double rho = rho_0 + beta * (salinity[i, j] - refSalinity);
                        newP_h[i, j] = j == gridPointsZ - 1 ? 0.0 : newP_h[i, j + 1] + rho * g * z_eta[i, j] * deta;

                        double drho_dx = beta * ds_dx;
                        double phys_z = z[i, j];
                        double baroclinicTerm = -(g / rho_0) * drho_dx * (estuaryDepth - phys_z);

                        // Advection in logical coordinates with toggle for third-order or first-order scheme
                        double advection_u_xi, advection_u_eta, advection_w_xi, advection_w_eta;
                        if (useThirdOrderAdvection)
                        {
                            advection_u_xi = u[i, j] > 0 ?
                                (i >= 2 ? (3 * u[i, j] - 4 * u[i - 1, j] + u[i - 2, j]) / (2 * dxi) : (u[i, j] - u[i - 1, j]) / dxi) :
                                (i < gridPointsX - 2 ? (u[i + 2, j] - 4 * u[i + 1, j] + 3 * u[i, j]) / (2 * dxi) : (u[i + 1, j] - u[i, j]) / dxi);
                            advection_u_eta = w[i, j] > 0 ?
                                (j >= 2 ? (3 * u[i, j] - 4 * u[i, j - 1] + u[i, j - 2]) / (2 * deta) : (u[i, j] - u[i, j - 1]) / deta) :
                                (j < gridPointsZ - 2 ? (u[i, j + 2] - 4 * u[i, j + 1] + 3 * u[i, j]) / (2 * deta) : (u[i, j + 1] - u[i, j]) / deta);
                            advection_w_xi = u[i, j] > 0 ?
                                (i >= 2 ? (3 * w[i, j] - 4 * w[i - 1, j] + w[i - 2, j]) / (2 * dxi) : (w[i, j] - w[i - 1, j]) / dxi) :
                                (i < gridPointsX - 2 ? (w[i + 2, j] - 4 * w[i + 1, j] + 3 * w[i, j]) / (2 * dxi) : (w[i + 1, j] - w[i, j]) / dxi);
                            advection_w_eta = w[i, j] > 0 ?
                                (j >= 2 ? (3 * w[i, j] - 4 * w[i, j - 1] + w[i, j - 2]) / (2 * deta) : (w[i, j] - w[i, j - 1]) / deta) :
                                (j < gridPointsZ - 2 ? (w[i, j + 2] - 4 * w[i, j + 1] + 3 * w[i, j]) / (2 * deta) : (w[i, j + 1] - w[i, j]) / deta);
                        }
                        else
                        {
                            advection_u_xi = u[i, j] > 0 ? (u[i, j] - u[i - 1, j]) / dxi : (u[i + 1, j] - u[i, j]) / dxi;
                            advection_u_eta = w[i, j] > 0 ? (u[i, j] - u[i, j - 1]) / deta : (u[i, j + 1] - u[i, j]) / deta;
                            advection_w_xi = u[i, j] > 0 ? (w[i, j] - w[i - 1, j]) / dxi : (w[i + 1, j] - w[i, j]) / dxi;
                            advection_w_eta = w[i, j] > 0 ? (w[i, j] - w[i, j - 1]) / deta : (w[i, j + 1] - w[i, j]) / deta;
                        }
                        double advection_u = (advection_u_xi * z_eta[i, j] - advection_u_eta * z_xi[i, j]) / J[i, j];
                        double advection_w = (advection_w_xi * z_eta[i, j] - advection_w_eta * z_xi[i, j]) / J[i, j];

                        // Diffusion in logical coordinates
                        double diffusion_u = newK_x[i, j] * ((u[i + 1, j] - 2 * u[i, j] + u[i - 1, j]) / (dxi * dxi) * (z_eta[i, j] * z_eta[i, j] + x_eta[i, j] * x_eta[i, j]) / J[i, j] +
                                            (u[i, j + 1] - 2 * u[i, j] + u[i, j - 1]) / (deta * deta) * (z_xi[i, j] * z_xi[i, j] + x_xi[i, j] * x_xi[i, j]) / J[i, j]) +
                                            newK_z[i, j] * ((u[i + 1, j] - 2 * u[i, j] + u[i - 1, j]) / (dxi * dxi) * (z_eta[i, j] * z_eta[i, j] + x_eta[i, j] * x_eta[i, j]) / J[i, j] +
                                            (u[i, j + 1] - 2 * u[i, j] + u[i, j - 1]) / (deta * deta) * (z_xi[i, j] * z_xi[i, j] + x_xi[i, j] * x_xi[i, j]) / J[i, j]);
                        double diffusion_w = newK_x[i, j] * ((w[i + 1, j] - 2 * w[i, j] + w[i - 1, j]) / (dxi * dxi) * (z_eta[i, j] * z_eta[i, j] + x_eta[i, j] * x_eta[i, j]) / J[i, j] +
                                            (w[i, j + 1] - 2 * w[i, j] + w[i, j - 1]) / (deta * deta) * (z_xi[i, j] * z_xi[i, j] + x_xi[i, j] * x_xi[i, j]) / J[i, j]) +
                                            newK_z[i, j] * ((w[i + 1, j] - 2 * w[i, j] + w[i - 1, j]) / (dxi * dxi) * (z_eta[i, j] * z_eta[i, j] + x_eta[i, j] * x_eta[i, j]) / J[i, j] +
                                            (w[i, j + 1] - 2 * w[i, j] + w[i, j - 1]) / (deta * deta) * (z_xi[i, j] * z_xi[i, j] + x_xi[i, j] * x_xi[i, j]) / J[i, j]);

                        uStar[i, j] = Math.Max(-2.0, Math.Min(2.0, u[i, j] + dt * (-advection_u + diffusion_u + baroclinicTerm)));
                        wStar[i, j] = Math.Max(-0.2, Math.Min(0.2, w[i, j] + dt * (-advection_w + diffusion_w + buoyancyValue)));
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
                            double divU = ((uStar[i + 1, j] - uStar[i - 1, j]) / (2 * dxi) * z_eta[i, j] - (uStar[i, j + 1] - uStar[i, j - 1]) / (2 * deta) * z_xi[i, j]) / J[i, j] +
                                          ((wStar[i + 1, j] - wStar[i - 1, j]) / (2 * dxi) * x_eta[i, j] - (wStar[i, j + 1] - wStar[i, j - 1]) / (2 * deta) * x_xi[i, j]) / J[i, j];
                            double laplacianP = ((p_nh[i + 1, j] - 2 * p_nh[i, j] + p_nh[i - 1, j]) / (dxi * dxi) * (z_eta[i, j] * z_eta[i, j] + x_eta[i, j] * x_eta[i, j]) +
                                                (p_nh[i, j + 1] - 2 * p_nh[i, j] + p_nh[i, j - 1]) / (deta * deta) * (z_xi[i, j] * z_xi[i, j] + x_xi[i, j] * x_xi[i, j])) / J[i, j];
                            newP_nh[i, j] = p_nh[i, j] + pressureRelaxation * (laplacianP - rho_0 * divU / dt);
                            newP_nh[i, j] = Math.Max(-1e4, Math.Min(1e4, newP_nh[i, j]));
                        }
                    }
                    p_nh = newP_nh;
                }

                // Step 3: Correct velocity with total pressure
                for (int i = 1; i < gridPointsX - 1; i++)
                {
                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        double dp_dxi = (p_h[i + 1, j] + p_nh[i + 1, j] - p_h[i - 1, j] - p_nh[i - 1, j]) / (2 * dxi);
                        double dp_deta = (p_h[i, j + 1] + p_nh[i, j + 1] - p_h[i, j - 1] - p_nh[i, j - 1]) / (2 * deta);
                        double dp_dx = (dp_dxi * z_eta[i, j] - dp_deta * z_xi[i, j]) / J[i, j];
                        double dp_dz = (-dp_dxi * x_eta[i, j] + dp_deta * x_xi[i, j]) / J[i, j];
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
                        double du_dxi = (u[i + 1, j] - u[i - 1, j]) / (2 * dxi);
                        double du_deta = (u[i, j + 1] - u[i, j - 1]) / (2 * deta);
                        double dw_dxi = (w[i + 1, j] - w[i - 1, j]) / (2 * dxi);
                        double dw_deta = (w[i, j + 1] - w[i, j - 1]) / (2 * deta);
                        double du_dx = (du_dxi * z_eta[i, j] - du_deta * z_xi[i, j]) / J[i, j];
                        double du_dz = (-du_dxi * x_eta[i, j] + du_deta * x_xi[i, j]) / J[i, j];
                        double dw_dx = (dw_dxi * z_eta[i, j] - dw_deta * z_xi[i, j]) / J[i, j];
                        double dw_dz = (-dw_dxi * x_eta[i, j] + dw_deta * x_xi[i, j]) / J[i, j];
                        double production = Math.Min(0.1, (newK_x[i, j] + newK_z[i, j]) / 2.0 *
                            (2 * du_dx * du_dx + 2 * dw_dz * dw_dz + (du_dz + dw_dx) * (du_dz + dw_dx)));

                        double salinityPerturbation = salinity[i, j] - refSalinity;
                        double buoyancyValue = Math.Max(-5e-2, Math.Min(5e-2, -g * beta * salinityPerturbation / rho_0));
                        newBuoyancy[i, j] = buoyancyValue;

                        double diffusion_k = (newK_x[i, j] / sigma_k) * ((k[i + 1, j] - 2 * k[i, j] + k[i - 1, j]) / (dxi * dxi) * (z_eta[i, j] * z_eta[i, j] + x_eta[i, j] * x_eta[i, j]) +
                            (k[i, j + 1] - 2 * k[i, j] + k[i, j - 1]) / (deta * deta) * (z_xi[i, j] * z_xi[i, j] + x_xi[i, j] * x_xi[i, j])) / J[i, j];
                        newK[i, j] = Math.Max(minK, k[i, j] + dt * (production - epsilon[i, j] + diffusion_k - buoyancyValue));

                        double diffusion_epsilon = (newK_x[i, j] / sigma_epsilon) * ((epsilon[i + 1, j] - 2 * epsilon[i, j] + epsilon[i - 1, j]) / (dxi * dxi) * (z_eta[i, j] * z_eta[i, j] + x_eta[i, j] * x_eta[i, j]) +
                            (epsilon[i, j + 1] - 2 * epsilon[i, j] + epsilon[i, j - 1]) / (deta * deta) * (z_xi[i, j] * z_xi[i, j] + x_xi[i, j] * x_xi[i, j])) / J[i, j];
                        double epsilonTerm = C1_epsilon * production * epsilon[i, j] / Math.Max(k[i, j], minK) -
                                            C2_epsilon * epsilon[i, j] * epsilon[i, j] / Math.Max(k[i, j], minK);
                        newEpsilon[i, j] = Math.Max(minEpsilon, epsilon[i, j] + dt * (epsilonTerm + diffusion_epsilon));

                        double ds_dxi = (salinity[i + 1, j] - salinity[i - 1, j]) / (2 * dxi);
                        double ds_deta = (salinity[i, j + 1] - salinity[i, j - 1]) / (2 * deta);
                        double advection_s_xi, advection_s_eta;
                        if (useThirdOrderAdvection)
                        {
                            advection_s_xi = u[i, j] > 0 ?
                                (i >= 2 ? (3 * salinity[i, j] - 4 * salinity[i - 1, j] + salinity[i - 2, j]) / (2 * dxi) : (salinity[i, j] - salinity[i - 1, j]) / dxi) :
                                (i < gridPointsX - 2 ? (salinity[i + 2, j] - 4 * salinity[i + 1, j] + 3 * salinity[i, j]) / (2 * dxi) : (salinity[i + 1, j] - salinity[i, j]) / dxi);
                            advection_s_eta = w[i, j] > 0 ?
                                (j >= 2 ? (3 * salinity[i, j] - 4 * salinity[i, j - 1] + salinity[i, j - 2]) / (2 * deta) : (salinity[i, j] - salinity[i, j - 1]) / deta) :
                                (j < gridPointsZ - 2 ? (salinity[i, j + 2] - 4 * salinity[i, j + 1] + 3 * salinity[i, j]) / (2 * deta) : (salinity[i, j + 1] - salinity[i, j]) / deta);
                        }
                        else
                        {
                            advection_s_xi = u[i, j] > 0 ? (salinity[i, j] - salinity[i - 1, j]) / dxi : (salinity[i + 1, j] - salinity[i, j]) / dxi;
                            advection_s_eta = w[i, j] > 0 ? (salinity[i, j] - salinity[i, j - 1]) / deta : (salinity[i, j + 1] - salinity[i, j]) / deta;
                        }
                        double advection_s = (advection_s_xi * z_eta[i, j] - advection_s_eta * z_xi[i, j]) / J[i, j];
                        double diffusion_s = newK_x[i, j] * ((salinity[i + 1, j] - 2 * salinity[i, j] + salinity[i - 1, j]) / (dxi * dxi) * (z_eta[i, j] * z_eta[i, j] + x_eta[i, j] * x_eta[i, j]) +
                            (salinity[i, j + 1] - 2 * salinity[i, j] + salinity[i, j - 1]) / (deta * deta) * (z_xi[i, j] * z_xi[i, j] + x_xi[i, j] * x_xi[i, j])) / J[i, j];
                        newSalinity[i, j] = Math.Max(0.0, Math.Min(maxSalinity, salinity[i, j] + dt * (diffusion_s - advection_s)));
                    }
                }

                // Boundary conditions
                for (int i = 0; i < gridPointsX; i++)
                {
                    double frictionVelocity = 0.05;
                    newK[i, 0] = newK[i, 1];
                    newK[i, gridPointsZ - 1] = newK[i, gridPointsZ - 2];
                    newEpsilon[i, 0] = Math.Pow(C_mu, 0.75) * Math.Pow(newK[i, 1], 1.5) / (vonKarman * z_eta[i, 1] * deta);
                    newEpsilon[i, gridPointsZ - 1] = newEpsilon[i, gridPointsZ - 2];
                    newSalinity[i, 0] = newSalinity[i, 1];
                    newSalinity[i, gridPointsZ - 1] = maxSalinity;
                    newBuoyancy[i, 0] = newBuoyancy[i, 1];
                    newBuoyancy[i, gridPointsZ - 1] = -g * beta * (newSalinity[i, gridPointsZ - 1] - refSalinity) / rho_0;
                    newP_h[i, 0] = newP_h[i, 1] + (rho_0 + beta * (newSalinity[i, 0] - refSalinity)) * g * z_eta[i, 1] * deta;
                    newP_h[i, gridPointsZ - 1] = 0.0;
                    p_nh[i, 0] = p_nh[i, 1];
                    p_nh[i, gridPointsZ - 1] = 0.0;
                    newK_x[i, 0] = newK_x[i, 1] * 0.1;
                    newK_x[i, gridPointsZ - 1] = newK_x[i, gridPointsZ - 2];
                    newK_z[i, 0] = newK_z[i, 1] * 0.01;
                    newK_z[i, gridPointsZ - 1] = newK_z[i, gridPointsZ - 2];
                    u[i, 0] = (frictionVelocity / vonKarman) * Math.Log(z[i, 1] / z0);
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
                    newP_h[0, j] = newP_h[1, j];
                    newP_h[gridPointsX - 1, j] = newP_h[gridPointsX - 2, j];
                    p_nh[0, j] = p_nh[1, j];
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
                    avgSalinity += salinity[i, j] * J[i, j];
                    avgK += k[i, j] * J[i, j];
                    avgEpsilon += epsilon[i, j] * J[i, j];
                    avgUSpeed += Math.Sqrt(u[i, j] * u[i, j] + w[i, j] * w[i, j]) * J[i, j];
                    avgK_x += K_x[i, j] * J[i, j];
                    avgK_z += K_z[i, j] * J[i, j];
                }
            double totalArea = 0.0;
            for (int i = 0; i < gridPointsX; i++)
                for (int j = 0; j < gridPointsZ; j++)
                    totalArea += J[i, j];
            avgSalinity /= totalArea;
            avgK /= totalArea;
            avgEpsilon /= totalArea;
            avgUSpeed /= totalArea;
            avgK_x /= totalArea;
            avgK_z /= totalArea;

            outputTextBox.AppendText($"Time: {time:F2}s | Avg Salinity: {avgSalinity:F2} PSU | Avg TKE: {avgK:E2} m²/s² | Avg Epsilon: {avgEpsilon:E2} m²/s³ | Avg Speed: {avgUSpeed:F4} m/s | Avg K_x: {avgK_x:F4} m²/s | Avg K_z: {avgK_z:F4} m²/s | dt: {dt:F6}s | Advection: {(useThirdOrderAdvection ? "Third-Order" : "First-Order")}\r\n");
        }

        private void VisualizationBox_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            int width = visualizationBox.Width;
            int height = visualizationBox.Height;
            double minSalinity = 0.0;
            double maxX = 0.0, maxZ = 0.0;
            for (int i = 0; i < gridPointsX; i++)
                for (int j = 0; j < gridPointsZ; j++)
                {
                    maxX = Math.Max(maxX, x[i, j]);
                    maxZ = Math.Max(maxZ, z[i, j]);
                }

            if (currentSlice == "Surface" || currentSlice == "Bottom")
            {
                int jIndex = currentSlice == "Surface" ? gridPointsZ - 1 : 0;
                for (int i = 0; i < gridPointsX - 1; i++)
                {
                    int px = (int)(x[i, jIndex] / maxX * width);
                    int pxNext = (int)(x[i + 1, jIndex] / maxX * width);
                    for (int pxCur = px; pxCur < pxNext && pxCur < width; pxCur++)
                    {
                        for (int py = 0; py < height; py++)
                        {
                            float normalized = (float)((salinity[i, jIndex] - minSalinity) / (maxSalinity - minSalinity));
                            normalized = Math.Max(0.0f, Math.Min(1.0f, normalized));
                            Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                            g.FillRectangle(new SolidBrush(color), pxCur, py, 1, 1);
                        }
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
            else if (currentSlice == "Buoyancy")
            {
                int jIndex = gridPointsZ - 1;
                for (int i = 0; i < gridPointsX - 1; i++)
                {
                    int px = (int)(x[i, jIndex] / maxX * width);
                    int pxNext = (int)(x[i + 1, jIndex] / maxX * width);
                    for (int pxCur = px; pxCur < pxNext && pxCur < width; pxCur++)
                    {
                        for (int py = 0; py < height; py++)
                        {
                            float normalized = (float)((buoyancy[i, jIndex] - minBuoyancy) / (maxBuoyancy - minBuoyancy));
                            normalized = Math.Max(0.0f, Math.Min(1.0f, normalized));
                            Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                            g.FillRectangle(new SolidBrush(color), pxCur, py, 1, 1);
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
                double maxSpeed = 0.0;
                for (int i = 0; i < gridPointsX; i++)
                    for (int j = 0; j < gridPointsZ; j++)
                        maxSpeed = Math.Max(maxSpeed, Math.Sqrt(u[i, j] * u[i, j] + w[i, j] * w[i, j]));
                for (int i = 0; i < gridPointsX; i += 5)
                {
                    for (int j = 0; j < gridPointsZ; j += 5)
                    {
                        int px = (int)(x[i, j] / maxX * width);
                        int py = height - (int)(z[i, j] / maxZ * height);
                        float speed = (float)Math.Sqrt(u[i, j] * u[i, j] + w[i, j] * w[i, j]);
                        float normalized = (float)(speed / (maxSpeed + 1e-6));
                        Color color = Color.FromArgb((int)(255 * normalized), 0, (int)(255 * (1 - normalized)));
                        g.FillRectangle(new SolidBrush(color), px, py, width / gridPointsX, height / gridPointsZ);
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
            else if (currentSlice == "Vertical Mixing")
            {
                double minK_z = minMixing;
                double maxK_z = maxMixing;
                for (int i = 0; i < gridPointsX; i++)
                {
                    for (int j = 0; j < gridPointsZ; j++)
                    {
                        int px = (int)(x[i, j] / maxX * width);
                        int py = height - (int)(z[i, j] / maxZ * height);
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
            else if (currentSlice == "Advection")
            {
                double maxAdvection = 0.0;
                double[,] advection_u = new double[gridPointsX, gridPointsZ];
                for (int i = 1; i < gridPointsX - 1; i++)
                {
                    for (int j = 1; j < gridPointsZ - 1; j++)
                    {
                        double advection_u_xi = useThirdOrderAdvection ?
                            (u[i, j] > 0 ? (i >= 2 ? (3 * u[i, j] - 4 * u[i - 1, j] + u[i - 2, j]) / (2 * dxi) : (u[i, j] - u[i - 1, j]) / dxi) :
                                           (i < gridPointsX - 2 ? (u[i + 2, j] - 4 * u[i + 1, j] + 3 * u[i, j]) / (2 * dxi) : (u[i + 1, j] - u[i, j]) / dxi)) :
                            (u[i, j] > 0 ? (u[i, j] - u[i - 1, j]) / dxi : (u[i + 1, j] - u[i, j]) / dxi);
                        double advection_u_eta = useThirdOrderAdvection ?
                            (w[i, j] > 0 ? (j >= 2 ? (3 * u[i, j] - 4 * u[i, j - 1] + u[i, j - 2]) / (2 * deta) : (u[i, j] - u[i, j - 1]) / deta) :
                                           (j < gridPointsZ - 2 ? (u[i, j + 2] - 4 * u[i, j + 1] + 3 * u[i, j]) / (2 * deta) : (u[i, j + 1] - u[i, j]) / deta)) :
                            (w[i, j] > 0 ? (u[i, j] - u[i, j - 1]) / deta : (u[i, j + 1] - u[i, j]) / deta);
                        advection_u[i, j] = Math.Abs((advection_u_xi * z_eta[i, j] - advection_u_eta * z_xi[i, j]) / J[i, j]);
                        maxAdvection = Math.Max(maxAdvection, advection_u[i, j]);
                    }
                }
                for (int i = 0; i < gridPointsX; i++)
                {
                    for (int j = 0; j < gridPointsZ; j++)
                    {
                        int px = (int)(x[i, j] / maxX * width);
                        int py = height - (int)(z[i, j] / maxZ * height);
                        float normalized = (float)(advection_u[i, j] / (maxAdvection + 1e-6));
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
                g.DrawString($"{maxAdvection:F6} m/s²", verdanaFont, Brushes.Black, width - 40, 0);
                g.DrawString("0 m/s²", verdanaFont, Brushes.Black, width - 40, height - 20);
            }
            else // Longitudinal
            {
                for (int i = 0; i < gridPointsX; i++)
                {
                    for (int j = 0; j < gridPointsZ; j++)
                    {
                        int px = (int)(x[i, j] / maxX * width);
                        int py = height - (int)(z[i, j] / maxZ * height);
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

        /**
        private void UseThirdOrderAdvectionCheckBox_CheckedChanged(object sender, EventArgs e)
        {
            useThirdOrderAdvection = useThirdOrderAdvectionCheckBox.Checked;
            if (isRunning)
            {
                simulationTimer.Stop();
                isRunning = false;
                InitializeSimulation();
                outputTextBox.AppendText($"Advection scheme changed to {(useThirdOrderAdvection ? "third-order" : "first-order")}. Simulation reset.\r\n");
                visualizationBox.Invalidate();
                UpdateButtonStates();
            }
            else
            {
                InitializeSimulation();
                outputTextBox.AppendText($"Advection scheme set to {(useThirdOrderAdvection ? "third-order" : "first-order")}.\r\n");
                visualizationBox.Invalidate();
            }
        }
        **/

        private void UseThirdOrderAdvectionCheckBox_CheckedChanged(object sender, EventArgs e)
        {
            useThirdOrderAdvection = useThirdOrderAdvectionCheckBox.Checked;
            simulationTimer.Stop();
            isRunning = false;
            InitializeSimulation();
            outputTextBox.AppendText($"Advection scheme changed to {(useThirdOrderAdvection ? "third-order" : "first-order")}. Simulation reset.\r\n");
            visualizationBox.Invalidate();
            UpdateButtonStates();
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

                deta = estuaryDepth / (gridPointsZ - 1);
                dt = ComputeTimeStep();
                simulationTimer.Start();
                isRunning = true;
                outputTextBox.AppendText($"Turbulence and Mixing simulation started with {(useThirdOrderAdvection ? "third-order" : "first-order")} advection.\r\n");
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
            outputTextBox.AppendText($"Simulation reset with {(useThirdOrderAdvection ? "third-order" : "first-order")} advection.\r\n");
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