using System;
using System.Numerics;
using System.Windows.Forms;
using System.Drawing;
using System.Collections.Generic;
using System.Linq;

namespace EstuarineCirculationModeling
{
    public class SpectralAnalyzer : Form
    {
        private static SpectralAnalyzer instance;
        private ComboBox variableComboBox;
        private TextBox durationTextBox;
        private TextBox samplingRateTextBox;
        private Button analyzeButton;
        private Button pauseButton;
        private Button resetButton;
        private Button clearButton;
        private Panel visualizationPanel;
        private TextBox outputTextBox;
        private Timer simulationTimer;
        private List<double> timeSeriesData;
        private double currentTime;
        private double duration;
        private double samplingRate;
        private bool isAnalyzing;

        // Custom class to hold spectral data
        private class SpectralData
        {
            public double[] PowerSpectrum { get; set; }
            public double[] Frequencies { get; set; }
        }

        public SpectralAnalyzer()
        {
            InitializeComponents();
            timeSeriesData = new List<double>();
            currentTime = 0;
            duration = 86400; // Default 1 day
            samplingRate = 3600; // Default 1 sample per hour
            isAnalyzing = false;
            InitializeSimulationTimer();
            UpdateButtonStates();
        }

        public static void ShowSpectralAnalyzerWindow()
        {
            if (instance == null || instance.IsDisposed)
            {
                instance = new SpectralAnalyzer();
            }
            instance.Show();
            instance.BringToFront();
        }

        private void InitializeComponents()
        {
            this.Text = "Spectral Analysis Tool";
            this.Size = new Size(800, 600);
            this.FormBorderStyle = FormBorderStyle.FixedDialog;
            this.MaximizeBox = false;
            this.Font = new Font("Consolas", 9F);

            // Control Panel
            Panel controlPanel = new Panel
            {
                Location = new Point(10, 10),
                Size = new Size(250, 500),
                BorderStyle = BorderStyle.FixedSingle,
                AutoScroll = true
            };

            // Variable Selection
            Label variableLabel = new Label
            {
                Text = "Select Variable:",
                Location = new Point(10, 10),
                AutoSize = true
            };

            variableComboBox = new ComboBox
            {
                Location = new Point(10, 30),
                Size = new Size(200, 22),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            variableComboBox.Items.AddRange(new string[] { "Richardson Number", "Velocity", "Water Level" });
            variableComboBox.SelectedIndex = 0;

            // Duration
            Label durationLabel = new Label
            {
                Text = "Duration (s):",
                Location = new Point(10, 60),
                AutoSize = true
            };

            durationTextBox = new TextBox
            {
                Location = new Point(10, 80),
                Size = new Size(200, 22),
                Text = "86400"
            };

            // Sampling Rate
            Label samplingRateLabel = new Label
            {
                Text = "Sampling Rate (s):",
                Location = new Point(10, 110),
                AutoSize = true
            };

            samplingRateTextBox = new TextBox
            {
                Location = new Point(10, 130),
                Size = new Size(200, 22),
                Text = "3600"
            };

            // Buttons
            analyzeButton = new Button
            {
                Text = "Analyze",
                Location = new Point(10, 160),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            analyzeButton.FlatAppearance.BorderSize = 1;
            analyzeButton.FlatAppearance.BorderColor = Color.Black;
            analyzeButton.Click += AnalyzeButton_Click;

            pauseButton = new Button
            {
                Text = "Pause",
                Location = new Point(120, 160),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            pauseButton.FlatAppearance.BorderSize = 1;
            pauseButton.FlatAppearance.BorderColor = Color.Black;
            pauseButton.Click += PauseButton_Click;

            resetButton = new Button
            {
                Text = "Reset",
                Location = new Point(10, 195),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            resetButton.FlatAppearance.BorderSize = 1;
            resetButton.FlatAppearance.BorderColor = Color.Black;
            resetButton.Click += ResetButton_Click;

            clearButton = new Button
            {
                Text = "Clear",
                Location = new Point(120, 195),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            clearButton.FlatAppearance.BorderSize = 1;
            clearButton.FlatAppearance.BorderColor = Color.Black;
            clearButton.Click += ClearButton_Click;

            // Visualization Panel
            visualizationPanel = new Panel
            {
                Location = new Point(270, 10),
                Size = new Size(500, 350),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            visualizationPanel.Paint += VisualizationPanel_Paint;

            // Output TextBox
            outputTextBox = new TextBox
            {
                Location = new Point(270, 370),
                Size = new Size(500, 150),
                Multiline = true,
                ReadOnly = true,
                ScrollBars = ScrollBars.Vertical,
                Font = new Font("Consolas", 9F)
            };

            // Add controls to control panel
            controlPanel.Controls.AddRange(new Control[] { variableLabel, variableComboBox, durationLabel, durationTextBox, samplingRateLabel, samplingRateTextBox, analyzeButton, pauseButton, resetButton, clearButton });

            // Add controls to form
            this.Controls.AddRange(new Control[] { controlPanel, visualizationPanel, outputTextBox });
        }

        private void InitializeSimulationTimer()
        {
            simulationTimer = new Timer();
            simulationTimer.Interval = 100; // Update every 100ms
            simulationTimer.Tick += (s, e) => UpdateSimulation();
        }

        private void UpdateButtonStates()
        {
            analyzeButton.Enabled = !isAnalyzing;
            pauseButton.Enabled = isAnalyzing;
            resetButton.Enabled = timeSeriesData.Count > 0;
            clearButton.Enabled = timeSeriesData.Count > 0 || visualizationPanel.Tag != null;
        }

        private void AnalyzeButton_Click(object sender, EventArgs e)
        {
            try
            {
                duration = double.Parse(durationTextBox.Text);
                samplingRate = double.Parse(samplingRateTextBox.Text);

                if (duration <= 0 || samplingRate <= 0)
                {
                    MessageBox.Show("Duration and sampling rate must be positive.", "Input Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                    return;
                }

                if (!isAnalyzing)
                {
                    timeSeriesData.Clear();
                    currentTime = 0;
                }
                isAnalyzing = true;
                simulationTimer.Start();
                outputTextBox.AppendText(isAnalyzing && currentTime == 0 ? "Starting spectral analysis simulation...\r\n" : "Resuming spectral analysis simulation...\r\n");
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
            isAnalyzing = false;
            outputTextBox.AppendText("Spectral analysis paused.\r\n");
            UpdateButtonStates();
        }

        private void ResetButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isAnalyzing = false;
            timeSeriesData.Clear();
            currentTime = 0;
            visualizationPanel.Tag = null;
            visualizationPanel.Invalidate();
            outputTextBox.AppendText("Spectral analysis reset.\r\n");
            UpdateButtonStates();
        }

        private void ClearButton_Click(object sender, EventArgs e)
        {
            simulationTimer.Stop();
            isAnalyzing = false;
            timeSeriesData.Clear();
            currentTime = 0;
            outputTextBox.Clear();
            visualizationPanel.Tag = null;
            visualizationPanel.Invalidate();
            outputTextBox.AppendText("Spectral analysis cleared.\r\n");
            UpdateButtonStates();
        }

        private void UpdateSimulation()
        {
            if (currentTime >= duration)
            {
                simulationTimer.Stop();
                isAnalyzing = false;
                PerformSpectralAnalysis();
                UpdateButtonStates();
                return;
            }

            // Simulate time series data based on selected variable
            string selectedVariable = variableComboBox.SelectedItem.ToString();
            double value = GenerateSyntheticData(selectedVariable, currentTime);
            timeSeriesData.Add(value);
            currentTime += samplingRate;
            visualizationPanel.Invalidate();
            UpdateButtonStates();
        }

        private double GenerateSyntheticData(string variable, double time)
        {
            // Simulate synthetic data for demonstration
            double tidalPeriod = 43200; // 12-hour tidal cycle
            double inertialPeriod = 17 * 3600; // Approx inertial period at 45° latitude
            double value = 0;

            if (variable == "Richardson Number")
            {
                // Simulate Ri with tidal and inertial components
                value = 0.5 + 0.3 * Math.Sin(2 * Math.PI * time / tidalPeriod) + 0.2 * Math.Sin(2 * Math.PI * time / inertialPeriod);
            }
            else if (variable == "Velocity")
            {
                // Simulate velocity with tidal and higher-frequency noise
                value = 0.1 + 0.05 * Math.Sin(2 * Math.PI * time / tidalPeriod) + 0.02 * Math.Sin(2 * Math.PI * time / (tidalPeriod / 2)) + new Random().NextDouble() * 0.01;
            }
            else if (variable == "Water Level")
            {
                // Simulate water level with tidal component
                value = 1.0 * Math.Sin(2 * Math.PI * time / tidalPeriod);
            }

            return value;
        }

        private void PerformSpectralAnalysis()
        {
            if (timeSeriesData.Count < 2)
            {
                outputTextBox.AppendText("Insufficient data for spectral analysis.\r\n");
                return;
            }

            // Perform FFT
            double[] data = timeSeriesData.ToArray();
            Complex[] complexData = new Complex[data.Length];
            for (int i = 0; i < data.Length; i++)
                complexData[i] = new Complex(data[i], 0);

            FourierTransform(complexData);

            // Calculate power spectrum
            double[] powerSpectrum = new double[complexData.Length / 2];
            double[] frequencies = new double[complexData.Length / 2];
            double deltaF = 1.0 / (samplingRate * complexData.Length);

            for (int i = 0; i < complexData.Length / 2; i++)
            {
                powerSpectrum[i] = complexData[i].Magnitude * complexData[i].Magnitude / complexData.Length;
                frequencies[i] = i * deltaF;
            }

            // Identify dominant frequencies
            List<(double Frequency, double Power)> dominantFrequencies = new List<(double, double)>();
            for (int i = 1; i < powerSpectrum.Length; i++)
            {
                if (powerSpectrum[i] > powerSpectrum.Max() * 0.1) // Threshold at 10% of max power
                {
                    dominantFrequencies.Add((frequencies[i], powerSpectrum[i]));
                }
            }

            // Output results
            outputTextBox.AppendText($"Spectral Analysis Results for {variableComboBox.SelectedItem}:\r\n");
            foreach (var (freq, power) in dominantFrequencies)
            {
                double period = 1.0 / freq;
                outputTextBox.AppendText($"Frequency: {freq:F6} Hz, Period: {period:F2} s, Power: {power:F4}\r\n");
                // Identify known physical periods
                if (Math.Abs(period - 43200) < 1000)
                    outputTextBox.AppendText("  -> Likely semi-diurnal tide (M2)\r\n");
                else if (Math.Abs(period - 17 * 3600) < 3600)
                    outputTextBox.AppendText("  -> Likely inertial oscillation\r\n");
            }

            visualizationPanel.Tag = new SpectralData { PowerSpectrum = powerSpectrum, Frequencies = frequencies };
            visualizationPanel.Invalidate();
        }

        private void VisualizationPanel_Paint(object sender, PaintEventArgs e)
        {
            Graphics g = e.Graphics;
            g.Clear(Color.White);
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

            const int margin = 50;
            float plotWidth = visualizationPanel.Width - 2 * margin;
            float plotHeight = visualizationPanel.Height - 2 * margin - 30; // Extra space for title
            float xOrigin = margin;
            float yOrigin = visualizationPanel.Height - margin;

            if (visualizationPanel.Tag is SpectralData spectralData)
            {
                double[] powerSpectrum = spectralData.PowerSpectrum;
                double[] frequencies = spectralData.Frequencies;
                double maxPower = powerSpectrum.Length > 0 ? powerSpectrum.Max() : 1.0;
                double maxFreq = frequencies.Length > 0 ? frequencies[frequencies.Length - 1] : 1.0;

                // Ensure maxPower is not too small to avoid division issues
                if (maxPower < 1E-10) maxPower = 1.0;

                // Draw title
                g.DrawString($"Power Spectrum of {variableComboBox.SelectedItem}", new Font("Consolas", 10F, FontStyle.Bold), Brushes.Black, margin, 10);

                // Draw axes
                g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin + plotWidth, yOrigin); // X-axis
                g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin, yOrigin - plotHeight); // Y-axis

                // Draw axis labels
                g.DrawString("Frequency (Hz)", new Font("Consolas", 9F), Brushes.Black, xOrigin + plotWidth / 2 - 30, yOrigin + 20);
                g.DrawString("Power", new Font("Consolas", 9F), Brushes.Black, xOrigin - 40, yOrigin - plotHeight / 2 - 10, new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical });

                // Draw grid and ticks
                DrawGridAndTicks(g, xOrigin, yOrigin, plotWidth, plotHeight, maxFreq, maxPower, true, margin);

                // Plot power spectrum
                for (int i = 1; i < powerSpectrum.Length; i++)
                {
                    float x1 = xOrigin + (float)(frequencies[i - 1] / maxFreq * plotWidth);
                    float y1 = yOrigin - (float)(powerSpectrum[i - 1] / maxPower * plotHeight);
                    float x2 = xOrigin + (float)(frequencies[i] / maxFreq * plotWidth);
                    float y2 = yOrigin - (float)(powerSpectrum[i] / maxPower * plotHeight);
                    // Clamp y1 and y2 to prevent overflow
                    y1 = Math.Max(margin, Math.Min(yOrigin, y1));
                    y2 = Math.Max(margin, Math.Min(yOrigin, y2));
                    g.DrawLine(Pens.Blue, x1, y1, x2, y2);
                }

                // Draw legend
                g.DrawLine(Pens.Blue, xOrigin + plotWidth - 60, 30, xOrigin + plotWidth - 40, 30);
                g.DrawString("Power Spectrum", new Font("Consolas", 8F), Brushes.Black, xOrigin + plotWidth - 100, 25);
            }
            else if (isAnalyzing)
            {
                double maxValue = timeSeriesData.Count > 0 ? Math.Max(Math.Abs(timeSeriesData.Min()), timeSeriesData.Max()) : 1.0;
                double maxTime = duration;

                // Ensure maxValue is not too small to avoid division issues
                if (maxValue < 1E-10) maxValue = 1.0;

                // Draw title
                g.DrawString($"Time Series: {variableComboBox.SelectedItem}", new Font("Consolas", 10F, FontStyle.Bold), Brushes.Black, margin, 10);

                // Draw axes
                g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin + plotWidth, yOrigin); // X-axis
                g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin, yOrigin - plotHeight); // Y-axis

                // Draw axis labels
                g.DrawString("Time (s)", new Font("Consolas", 9F), Brushes.Black, xOrigin + plotWidth / 2 - 20, yOrigin + 20);
                g.DrawString(variableComboBox.SelectedItem.ToString(), new Font("Consolas", 9F), Brushes.Black, xOrigin - 40, yOrigin - plotHeight / 2 - 10, new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical });

                // Draw grid and ticks
                DrawGridAndTicks(g, xOrigin, yOrigin, plotWidth, plotHeight, maxTime, maxValue, false, margin);

                // Plot time series
                for (int i = 1; i < timeSeriesData.Count; i++)
                {
                    float x1 = xOrigin + (float)((i - 1) * samplingRate / maxTime * plotWidth);
                    float y1 = yOrigin - (float)(timeSeriesData[i - 1] / maxValue * plotHeight / 2 + plotHeight / 2);
                    float x2 = xOrigin + (float)(i * samplingRate / maxTime * plotWidth);
                    float y2 = yOrigin - (float)(timeSeriesData[i] / maxValue * plotHeight / 2 + plotHeight / 2);
                    // Clamp y1 and y2 to prevent overflow
                    y1 = Math.Max(margin, Math.Min(yOrigin, y1));
                    y2 = Math.Max(margin, Math.Min(yOrigin, y2));
                    g.DrawLine(Pens.Blue, x1, y1, x2, y2);
                }

                // Draw legend
                g.DrawLine(Pens.Blue, xOrigin + plotWidth - 60, 30, xOrigin + plotWidth - 40, 30);
                g.DrawString("Time Series", new Font("Consolas", 8F), Brushes.Black, xOrigin + plotWidth - 100, 25);
            }
        }

        private void DrawGridAndTicks(Graphics g, float xOrigin, float yOrigin, float plotWidth, float plotHeight, double maxX, double maxY, bool isPowerSpectrum, int margin)
        {
            // Ensure maxX and maxY are not too small to avoid division issues
            if (maxX < 1E-10) maxX = 1.0;
            if (maxY < 1E-10) maxY = 1.0;

            // Calculate tick intervals
            double xTickInterval = CalculateTickInterval(maxX, 5); // Aim for ~5 ticks
            double yTickInterval = CalculateTickInterval(maxY, 5);

            // Draw vertical grid lines and x-axis ticks
            for (double x = 0; x <= maxX; x += xTickInterval)
            {
                float xPos = xOrigin + (float)(x / maxX * plotWidth);
                g.DrawLine(Pens.LightGray, xPos, yOrigin, xPos, yOrigin - plotHeight); // Grid line
                g.DrawLine(Pens.Black, xPos, yOrigin, xPos, yOrigin + 5); // Tick
                string label = isPowerSpectrum ? x.ToString("F4") : ((int)x).ToString();
                g.DrawString(label, new Font("Consolas", 8F), Brushes.Black, xPos - 15, yOrigin + 10);
            }

            // Draw horizontal grid lines and y-axis ticks
            double yStart = isPowerSpectrum ? 0 : -maxY;
            for (double y = yStart; y <= maxY; y += yTickInterval)
            {
                float yPos = yOrigin - (float)(y / maxY * plotHeight / (isPowerSpectrum ? 1 : 2) + (isPowerSpectrum ? 0 : plotHeight / 2));
                // Clamp yPos to prevent overflow
                yPos = Math.Max(margin, Math.Min(yOrigin, yPos));
                g.DrawLine(Pens.LightGray, xOrigin, yPos, xOrigin + plotWidth, yPos); // Grid line
                g.DrawLine(Pens.Black, xOrigin - 5, yPos, xOrigin, yPos); // Tick
                string label = y.ToString("F2");
                g.DrawString(label, new Font("Consolas", 8F), Brushes.Black, xOrigin - 40, yPos - 5);
            }
        }

        private double CalculateTickInterval(double range, int targetTicks)
        {
            // Calculate a reasonable tick interval based on the range
            double roughInterval = range / targetTicks;
            double[] niceNumbers = { 1.0, 2.0, 5.0, 10.0 };
            double magnitude = Math.Pow(10, Math.Floor(Math.Log10(roughInterval)));
            double normalized = roughInterval / magnitude;
            double niceInterval = niceNumbers[0];
            for (int i = 1; i < niceNumbers.Length; i++)
            {
                if (normalized <= niceNumbers[i])
                {
                    niceInterval = niceNumbers[i];
                    break;
                }
            }
            return niceInterval * magnitude;
        }

        private void FourierTransform(Complex[] data)
        {
            int N = data.Length;
            if (N <= 1) return;

            // Divide
            Complex[] even = new Complex[N / 2];
            Complex[] odd = new Complex[N / 2];
            for (int i = 0; i < N / 2; i++)
            {
                even[i] = data[2 * i];
                odd[i] = data[2 * i + 1];
            }

            // Conquer
            FourierTransform(even);
            FourierTransform(odd);

            // Combine
            for (int k = 0; k < N / 2; k++)
            {
                Complex t = Complex.Exp(-2 * Math.PI * Complex.ImaginaryOne * k / N) * odd[k];
                data[k] = even[k] + t;
                data[k + N / 2] = even[k] - t;
            }
        }
    }
}