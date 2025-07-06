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
        private Stratification stratification; // Add Stratification instance

        public Form1()
        {
            InitializeComponent();
            model = new EstuarineModel();
            stratification = new Stratification(model.EstuaryLength, model.EstuaryDepth, 100); // Initialize with model parameters
            isSimulationRunning = false;
            InitializeSimulationTimer();
            UpdateButtonStates();
            this.Font = new Font("Tahoma", 8.25F);
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
                model.RiverInflow = double.Parse(riverInflowTextBox.Text);
                model.TidalAmplitude = double.Parse(tidalAmplitudeTextBox.Text);
                model.TidalPeriod = double.Parse(tidalPeriodTextBox.Text);
                model.SalinityOcean = double.Parse(oceanSalinityTextBox.Text);
                model.TemperatureOcean = double.Parse(oceanTemperatureTextBox.Text);
                model.EstuaryLength = double.Parse(estuaryLengthTextBox.Text);
                model.EstuaryDepth = double.Parse(estuaryDepthTextBox.Text);
                model.UseRANSSolver = useRANSSolverCheckBox.Checked;
                model.UseNonHydrostaticOverride = useNonHydrostaticCheckBox.Checked;
                // Reinitialize stratification with updated parameters
                stratification = new Stratification(model.EstuaryLength, model.EstuaryDepth, 100);
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
            stratificationButton.Enabled = model.UseRANSSolver; // Enable stratification button with RANS
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

        private void UpdateButtonStates()
        {
            startButton.Enabled = !isSimulationRunning;
            pauseButton.Enabled = isSimulationRunning;
            resetButton.Enabled = true;
            useNonHydrostaticCheckBox.Enabled = model.UseRANSSolver && !isSimulationRunning;
            stratificationButton.Enabled = model.UseRANSSolver && !isSimulationRunning; // Stratification button state
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
            VisualizationRenderer.Render(e.Graphics, visualizationPanel.Width, visualizationPanel.Height, model);
        }
    }
}