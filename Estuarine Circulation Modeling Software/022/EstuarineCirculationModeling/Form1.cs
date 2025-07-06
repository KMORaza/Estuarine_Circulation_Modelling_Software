using System;
using System.Windows.Forms;
using System.Drawing;

namespace EstuarineCirculationModeling
{
    public partial class Form1 : Form
    {
        private EstuarineModel model;
        private Timer simulationTimer;
        private bool isSimulationRunning;
        private Stratification stratification;
        private PassiveScalarTransportEq passiveScalarTransport;
        private TidalStrainSim tidalStrainSim;
        private TurbulenceAndMixing turbulenceAndMixing;
        private double riverScalarConcentration = 1.0; // Default river scalar concentration (kg/m³)

        public Form1()
        {
            InitializeComponent();
            model = new EstuarineModel();
            stratification = new Stratification(model.EstuaryLength, model.EstuaryDepth, 100);
            passiveScalarTransport = new PassiveScalarTransportEq(model.EstuaryLength, 100);
            tidalStrainSim = new TidalStrainSim(model.EstuaryLength, model.EstuaryDepth, 43200, 1.0);
            turbulenceAndMixing = new TurbulenceAndMixing();
            isSimulationRunning = false;
            InitializeSimulationTimer();
            UpdateButtonStates();
            this.Font = new Font("Verdana", 8.25F);
        }

        private void InitializeSimulationTimer()
        {
            simulationTimer = new Timer();
            simulationTimer.Interval = 100; // Update every 100ms
            simulationTimer.Tick += (s, e) => UpdateSimulation();
        }

        private void UpdateSimulation()
        {
            model.Update();
            // Apply stratification effects if RANS solver is used
            if (model.UseRANSSolver)
            {
                double[] velocityProfile = model.GetVelocityProfile();
                double[] eddyViscosity = new double[velocityProfile.Length];
                for (int i = 0; i < eddyViscosity.Length; i++)
                    eddyViscosity[i] = model.GetEddyViscosityAtPoint(i * model.EstuaryLength / velocityProfile.Length);
                stratification.ComputeStratification(velocityProfile, model.GetSalinityProfile(), model.GetTemperatureProfile(), eddyViscosity);
            }
            visualizationPanel.Invalidate();
            UpdateOutputConsole();
        }

        private void startButton_Click(object sender, EventArgs e)
        {
            try
            {
                double riverInflow = double.Parse(riverInflowTextBox.Text);
                double tidalAmplitude = double.Parse(tidalAmplitudeTextBox.Text);
                double tidalPeriod = double.Parse(tidalPeriodTextBox.Text);
                double salinityOcean = double.Parse(oceanSalinityTextBox.Text);
                double temperatureOcean = double.Parse(oceanTemperatureTextBox.Text);
                double estuaryLength = double.Parse(estuaryLengthTextBox.Text);
                double estuaryDepth = double.Parse(estuaryDepthTextBox.Text);

                // Clamp inputs to prevent invalid values
                model.RiverInflow = Math.Max(0.01, Math.Min(100.0, riverInflow));
                model.TidalAmplitude = Math.Max(0.1, Math.Min(10.0, tidalAmplitude));
                model.TidalPeriod = Math.Max(3600.0, Math.Min(86400.0, tidalPeriod));
                model.SalinityOcean = Math.Max(0.1, Math.Min(40.0, salinityOcean)); // Prevent zero
                model.TemperatureOcean = Math.Max(0.1, Math.Min(40.0, temperatureOcean)); // Prevent zero
                model.EstuaryLength = Math.Max(1000.0, Math.Min(100000.0, estuaryLength));
                model.EstuaryDepth = Math.Max(1.0, Math.Min(100.0, estuaryDepth));
                model.UseRANSSolver = useRANSSolverCheckBox.Checked;
                model.UseNonHydrostaticOverride = useNonHydrostaticCheckBox.Checked;

                // Reinitialize stratification, passive scalar transport, and tidal strain with updated parameters
                stratification = new Stratification(model.EstuaryLength, model.EstuaryDepth, 100);
                passiveScalarTransport = new PassiveScalarTransportEq(model.EstuaryLength, 100);
                tidalStrainSim = new TidalStrainSim(model.EstuaryLength, model.EstuaryDepth, model.TidalPeriod, model.TidalAmplitude);

                outputConsoleTextBox.AppendText($"Simulation started with {(model.UseRANSSolver ? "RANS" : "simple")} solver ({(model.IsNonHydrostatic() ? "Non-Hydrostatic" : "Hydrostatic")}).\r\n");
                model.Reset();
                simulationTimer.Start();
                isSimulationRunning = true;
                UpdateButtonStates();
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
            outputConsoleTextBox.AppendText("Simulation paused.\r\n");
            UpdateButtonStates();
        }

        private void resetButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isSimulationRunning = false;
            model.Reset();
            outputConsoleTextBox.Clear();
            outputConsoleTextBox.AppendText("Simulation reset.\r\n");
            visualizationPanel.Invalidate();
            UpdateButtonStates();
        }

        private void useRANSSolverCheckBox_CheckedChanged(object sender, EventArgs e)
        {
            model.UseRANSSolver = useRANSSolverCheckBox.Checked;
            useNonHydrostaticCheckBox.Enabled = model.UseRANSSolver;
            stratificationButton.Enabled = model.UseRANSSolver;
            largeEddyButton.Enabled = true; // LES button always enabled
            lbLesButton.Enabled = true; // LB-LES button always enabled
            tidalStrainButton.Enabled = true; // Tidal Strain button always enabled
            satButton.Enabled = true; // SAT button always enabled
            compForcingButton.Enabled = true; // Comp. Forcing button always enabled
            tamButton.Enabled = true; // TAM button always enabled
            outputConsoleTextBox.AppendText($"Switched to {(model.UseRANSSolver ? "RANS" : "simple")} solver.\r\n");
        }

        private void useNonHydrostaticCheckBox_CheckedChanged(object sender, EventArgs e)
        {
            model.UseNonHydrostaticOverride = useNonHydrostaticCheckBox.Checked;
            outputConsoleTextBox.AppendText($"Switched to {(model.IsNonHydrostatic() ? "non-hydrostatic" : "hydrostatic")} approximation.\r\n");
        }

        private void stratificationButton_Click(object sender, EventArgs e)
        {
            stratification.ShowStratificationWindow();
        }

        private void largeEddyButton_Click(object sender, EventArgs e)
        {
            LargeEddySim.ShowLargeEddyWindow();
        }

        private void lbLesButton_Click(object sender, EventArgs e)
        {
            LatticeBoltzmannLES.ShowLBLESWindow();
        }

        private void tidalStrainButton_Click(object sender, EventArgs e)
        {
            tidalStrainSim.ShowTidalStrainWindow();
        }

        private void satButton_Click(object sender, EventArgs e)
        {
            SpectralAnalyzer.ShowSpectralAnalyzerWindow();
        }

        private void compForcingButton_Click(object sender, EventArgs e)
        {
            CompForcingMechanism.ShowCompForcingWindow();
        }

        private void tamButton_Click(object sender, EventArgs e)
        {
            turbulenceAndMixing.ShowTurbulenceWindow();
        }

        private void UpdateButtonStates()
        {
            startButton.Enabled = !isSimulationRunning;
            pauseButton.Enabled = isSimulationRunning;
            resetButton.Enabled = true;
            useNonHydrostaticCheckBox.Enabled = model.UseRANSSolver && !isSimulationRunning;
            stratificationButton.Enabled = model.UseRANSSolver && !isSimulationRunning;
            largeEddyButton.Enabled = true; // LES button always enabled
            lbLesButton.Enabled = true; // LB-LES button always enabled
            tidalStrainButton.Enabled = true; // Tidal Strain button always enabled
            satButton.Enabled = true; // SAT button always enabled
            compForcingButton.Enabled = true; // Comp. Forcing button always enabled
            tamButton.Enabled = true; // TAM button always enabled
        }

        private void UpdateOutputConsole()
        {
            double waterLevel = model.TidalAmplitude * Math.Sin(2 * Math.PI * model.CurrentTime / model.TidalPeriod);
            double temperatureAtWedge = model.GetTemperatureAtPoint(model.SaltWedgePosition);
            string output = $"Time: {model.CurrentTime:F2}s | Salt Wedge Position: {model.SaltWedgePosition:F2}m | Max Salinity: {model.GetMaxSalinity():F2} PSU | Temperature: {temperatureAtWedge:F2} °C | Water Level: {waterLevel:F2}m";
            if (model.UseRANSSolver)
            {
                double[] velocityProfile = model.GetVelocityProfile();
                double avgVelocity = 0.0;
                for (int i = 0; i < velocityProfile.Length; i++)
                    avgVelocity += velocityProfile[i];
                avgVelocity /= velocityProfile.Length;
                double[] richardsonNumber = stratification.GetRichardsonNumber();
                double avgRi = 0.0;
                for (int i = 0; i < richardsonNumber.Length; i++)
                    avgRi += richardsonNumber[i];
                avgRi /= richardsonNumber.Length;
                output += $" | Avg Velocity: {avgVelocity:F4} m/s | Avg Ri: {avgRi:F4} | {(model.IsNonHydrostatic() ? "Non-Hydrostatic" : "Hydrostatic")}";
            }
            outputConsoleTextBox.AppendText(output + "\r\n");
        }

        private void visualizationPanel_Paint(object sender, PaintEventArgs e)
        {
            VisualizationRenderer.Render(e.Graphics, visualizationPanel.Width, visualizationPanel.Height, model, stratification, passiveScalarTransport, riverScalarConcentration);
        }
    }
}