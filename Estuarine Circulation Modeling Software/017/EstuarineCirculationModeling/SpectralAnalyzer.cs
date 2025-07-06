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
        private ComboBox scaleComboBox;
        private ComboBox windowComboBox;
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
        private List<MixingEvent> mixingEvents;
        private double currentTime;
        private double duration;
        private double samplingRate;
        private bool isAnalyzing;
        private string selectedWindowType;

        // Custom class to hold spectral data
        private class SpectralData
        {
            public double[] PowerSpectrum { get; set; }
            public double[] Frequencies { get; set; }
            public int SegmentLength { get; set; }
        }

        // Custom class to hold mixing event data
        private class MixingEvent
        {
            public double Time { get; set; }
            public double Value { get; set; }
        }

        public SpectralAnalyzer()
        {
            InitializeComponents();
            timeSeriesData = new List<double>();
            mixingEvents = new List<MixingEvent>();
            currentTime = 0;
            duration = 86400; // Default 1 day
            samplingRate = 3600; // Default 1 sample per hour
            isAnalyzing = false;
            selectedWindowType = "Hanning"; // Default window
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
            this.FormClosing += SpectralAnalyzer_FormClosing;

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

            // Scale Selection
            Label scaleLabel = new Label
            {
                Text = "Plot Scale:",
                Location = new Point(10, 60),
                AutoSize = true
            };

            scaleComboBox = new ComboBox
            {
                Location = new Point(10, 80),
                Size = new Size(200, 22),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            scaleComboBox.Items.AddRange(new string[] { "Linear", "Log-Log" });
            scaleComboBox.SelectedIndex = 0;
            scaleComboBox.SelectedIndexChanged += (s, e) => { if (!IsDisposed) visualizationPanel.Invalidate(); };

            // Window Selection
            Label windowLabel = new Label
            {
                Text = "Window Type:",
                Location = new Point(10, 110),
                AutoSize = true
            };

            windowComboBox = new ComboBox
            {
                Location = new Point(10, 130),
                Size = new Size(200, 22),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            windowComboBox.Items.AddRange(new string[] { "Hanning", "Hamming", "Blackman", "Rectangular" });
            windowComboBox.SelectedIndex = 0; // Default to Hanning
            windowComboBox.SelectedIndexChanged += (s, e) =>
            {
                if (!IsDisposed)
                {
                    selectedWindowType = windowComboBox.SelectedItem?.ToString() ?? "Hanning";
                    visualizationPanel.Invalidate(); // Refresh plot to reflect window change in title
                }
            };

            // Duration
            Label durationLabel = new Label
            {
                Text = "Duration (s):",
                Location = new Point(10, 160),
                AutoSize = true
            };

            durationTextBox = new TextBox
            {
                Location = new Point(10, 180),
                Size = new Size(200, 22),
                Text = "86400"
            };

            // Sampling Rate
            Label samplingRateLabel = new Label
            {
                Text = "Sampling Rate (s):",
                Location = new Point(10, 210),
                AutoSize = true
            };

            samplingRateTextBox = new TextBox
            {
                Location = new Point(10, 230),
                Size = new Size(200, 22),
                Text = "3600"
            };

            // Buttons
            analyzeButton = new Button
            {
                Text = "Analyze",
                Location = new Point(10, 260),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            analyzeButton.FlatAppearance.BorderSize = 1;
            analyzeButton.FlatAppearance.BorderColor = Color.Black;
            analyzeButton.Click += AnalyzeButton_Click;

            pauseButton = new Button
            {
                Text = "Pause",
                Location = new Point(120, 260),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            pauseButton.FlatAppearance.BorderSize = 1;
            pauseButton.FlatAppearance.BorderColor = Color.Black;
            pauseButton.Click += PauseButton_Click;

            resetButton = new Button
            {
                Text = "Reset",
                Location = new Point(10, 295),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            resetButton.FlatAppearance.BorderSize = 1;
            resetButton.FlatAppearance.BorderColor = Color.Black;
            resetButton.Click += ResetButton_Click;

            clearButton = new Button
            {
                Text = "Clear",
                Location = new Point(120, 295),
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
            controlPanel.Controls.AddRange(new Control[] { variableLabel, variableComboBox, scaleLabel, scaleComboBox, windowLabel, windowComboBox, durationLabel, durationTextBox, samplingRateLabel, samplingRateTextBox, analyzeButton, pauseButton, resetButton, clearButton });

            // Add controls to form
            this.Controls.AddRange(new Control[] { controlPanel, visualizationPanel, outputTextBox });
        }

        private void SpectralAnalyzer_FormClosing(object sender, FormClosingEventArgs e)
        {
            simulationTimer.Stop();
            isAnalyzing = false;
        }

        private void InitializeSimulationTimer()
        {
            simulationTimer = new Timer();
            simulationTimer.Interval = 100; // Update every 100ms
            simulationTimer.Tick += (s, e) => UpdateSimulation();
        }

        private void UpdateButtonStates()
        {
            if (IsDisposed) return;

            if (InvokeRequired)
            {
                Invoke((Action)UpdateButtonStates);
                return;
            }

            analyzeButton.Enabled = !isAnalyzing;
            pauseButton.Enabled = isAnalyzing;
            resetButton.Enabled = timeSeriesData.Count > 0;
            clearButton.Enabled = timeSeriesData.Count > 0 || visualizationPanel.Tag != null;
        }

        private void AnalyzeButton_Click(object sender, EventArgs e)
        {
            if (IsDisposed) return;

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
                    mixingEvents.Clear();
                    currentTime = 0;
                }
                isAnalyzing = true;
                simulationTimer.Start();
                if (InvokeRequired)
                {
                    Invoke((Action)(() => outputTextBox.AppendText(isAnalyzing && currentTime == 0 ? "Starting spectral analysis simulation...\r\n" : "Resuming spectral analysis simulation...\r\n")));
                }
                else
                {
                    outputTextBox.AppendText(isAnalyzing && currentTime == 0 ? "Starting spectral analysis simulation...\r\n" : "Resuming spectral analysis simulation...\r\n");
                }
                UpdateButtonStates();
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Error: {ex.Message}", "Input Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
            }
        }

        private void PauseButton_Click(object sender, EventArgs e)
        {
            if (IsDisposed) return;

            simulationTimer.Stop();
            isAnalyzing = false;
            if (InvokeRequired)
            {
                Invoke((Action)(() => outputTextBox.AppendText("Spectral analysis paused.\r\n")));
            }
            else
            {
                outputTextBox.AppendText("Spectral analysis paused.\r\n");
            }
            UpdateButtonStates();
        }

        private void ResetButton_Click(object sender, EventArgs e)
        {
            if (IsDisposed) return;

            simulationTimer.Stop();
            isAnalyzing = false;
            timeSeriesData.Clear();
            mixingEvents.Clear();
            currentTime = 0;
            visualizationPanel.Tag = null;
            if (!visualizationPanel.IsDisposed) visualizationPanel.Invalidate();
            if (InvokeRequired)
            {
                Invoke((Action)(() => outputTextBox.AppendText("Spectral analysis reset.\r\n")));
            }
            else
            {
                outputTextBox.AppendText("Spectral analysis reset.\r\n");
            }
            UpdateButtonStates();
        }

        private void ClearButton_Click(object sender, EventArgs e)
        {
            if (IsDisposed) return;

            simulationTimer.Stop();
            isAnalyzing = false;
            timeSeriesData.Clear();
            mixingEvents.Clear();
            currentTime = 0;
            if (InvokeRequired)
            {
                Invoke((Action)(() => outputTextBox.Clear()));
                Invoke((Action)(() => outputTextBox.AppendText("Spectral analysis cleared.\r\n")));
            }
            else
            {
                outputTextBox.Clear();
                outputTextBox.AppendText("Spectral analysis cleared.\r\n");
            }
            visualizationPanel.Tag = null;
            if (!visualizationPanel.IsDisposed) visualizationPanel.Invalidate();
            UpdateButtonStates();
        }

        private void UpdateSimulation()
        {
            if (IsDisposed) return;

            if (currentTime >= duration)
            {
                simulationTimer.Stop();
                isAnalyzing = false;
                PerformSpectralAnalysis();
                UpdateButtonStates();
                return;
            }

            // Simulate time series data based on selected variable
            string selectedVariable = variableComboBox.SelectedItem?.ToString() ?? "Unknown";
            double value = GenerateSyntheticData(selectedVariable, currentTime);
            timeSeriesData.Add(value);

            // Detect mixing events (high variability)
            if (timeSeriesData.Count >= 5) // Need enough points for standard deviation
            {
                int windowSize = Math.Min(5, timeSeriesData.Count);
                var window = timeSeriesData.Skip(timeSeriesData.Count - windowSize).Take(windowSize).ToList();
                double mean = window.Average();
                double variance = window.Sum(x => (x - mean) * (x - mean)) / windowSize;
                double stdDev = Math.Sqrt(variance);
                double threshold = selectedVariable == "Richardson Number" ? 0.2 : selectedVariable == "Velocity" ? 0.03 : 0.5;
                if (stdDev > threshold)
                {
                    mixingEvents.Add(new MixingEvent { Time = currentTime, Value = value });
                }
            }

            currentTime += samplingRate;
            if (!visualizationPanel.IsDisposed) visualizationPanel.Invalidate();
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
                // Simulate Ri with tidal, inertial components, and occasional spikes for mixing
                value = 0.5 + 0.3 * Math.Sin(2 * Math.PI * time / tidalPeriod) + 0.2 * Math.Sin(2 * Math.PI * time / inertialPeriod);
                if (new Random((int)time).NextDouble() < 0.02) value += 0.5; // Simulate mixing event
            }
            else if (variable == "Velocity")
            {
                // Simulate velocity with tidal, higher-frequency noise, and turbulence-like fluctuations
                value = 0.1 + 0.05 * Math.Sin(2 * Math.PI * time / tidalPeriod) + 0.02 * Math.Sin(2 * Math.PI * time / (tidalPeriod / 2)) + new Random((int)time).NextDouble() * 0.01;
            }
            else if (variable == "Water Level")
            {
                // Simulate water level with M2 and K1 tidal components
                value = 1.0 * Math.Sin(2 * Math.PI * time / 44712) + 0.5 * Math.Sin(2 * Math.PI * time / 86148); // M2 + K1
            }

            return value;
        }

        private void PerformSpectralAnalysis()
        {
            if (IsDisposed) return;

            if (timeSeriesData.Count < 2)
            {
                if (InvokeRequired)
                {
                    Invoke((Action)(() => outputTextBox.AppendText("Insufficient data for spectral analysis.\r\n")));
                }
                else
                {
                    outputTextBox.AppendText("Insufficient data for spectral analysis.\r\n");
                }
                return;
            }

            // Perform Welch's method for PSD
            double[] data = timeSeriesData.ToArray();
            int segmentLength = Math.Min(256, data.Length);
            int overlap = segmentLength / 2;
            int numSegments = (data.Length - overlap) / (segmentLength - overlap);
            if (numSegments < 1) numSegments = 1;

            double[] window = new double[segmentLength];
            double windowSum = 0;
            for (int i = 0; i < segmentLength; i++)
            {
                double n = 2.0 * Math.PI * i / (segmentLength - 1);
                if (selectedWindowType == "Hanning")
                {
                    window[i] = 0.5 * (1 - Math.Cos(n));
                }
                else if (selectedWindowType == "Hamming")
                {
                    window[i] = 0.54 - 0.46 * Math.Cos(n);
                }
                else if (selectedWindowType == "Blackman")
                {
                    window[i] = 0.42 - 0.5 * Math.Cos(n) + 0.08 * Math.Cos(2 * n);
                }
                else // Rectangular
                {
                    window[i] = 1.0;
                }
                windowSum += window[i] * window[i];
            }
            double windowPower = windowSum / segmentLength;

            int fftLength = segmentLength / 2;
            double[] psd = new double[fftLength];
            double[] frequencies = new double[fftLength];
            double deltaF = 1.0 / (samplingRate * segmentLength);

            for (int i = 0; i < fftLength; i++)
            {
                frequencies[i] = i * deltaF;
                psd[i] = 0;
            }

            for (int seg = 0; seg < numSegments; seg++)
            {
                int start = seg * (segmentLength - overlap);
                Complex[] segmentData = new Complex[segmentLength];
                for (int i = 0; i < segmentLength; i++)
                {
                    double windowedValue = (start + i < data.Length) ? data[start + i] * window[i] : 0;
                    segmentData[i] = new Complex(windowedValue, 0);
                }

                FourierTransform(segmentData);

                for (int i = 0; i < fftLength; i++)
                {
                    psd[i] += segmentData[i].Magnitude * segmentData[i].Magnitude;
                }
            }

            // Average and normalize PSD
            for (int i = 0; i < fftLength; i++)
            {
                psd[i] = psd[i] / (numSegments * windowPower * samplingRate);
            }

            // Identify dominant frequencies and tidal constituents
            List<(double Frequency, double Power, string Label)> dominantFrequencies = new List<(double, double, string)>();
            double m2Freq = 1.0 / 44712; // M2: ~12.42 hours
            double k1Freq = 1.0 / 86148; // K1: ~23.93 hours
            double freqTolerance = deltaF * 2; // Tolerance for matching

            for (int i = 1; i < psd.Length; i++)
            {
                string label = "";
                if (Math.Abs(frequencies[i] - m2Freq) < freqTolerance)
                    label = "M2 (Semi-diurnal)";
                else if (Math.Abs(frequencies[i] - k1Freq) < freqTolerance)
                    label = "K1 (Diurnal)";
                else if (psd[i] > psd.Max() * 0.1)
                    label = "Other";

                if (!string.IsNullOrEmpty(label))
                {
                    dominantFrequencies.Add((frequencies[i], psd[i], label));
                }
            }

            // Analyze turbulence (check for -5/3 spectral slope)
            bool turbulenceDetected = false;
            if (frequencies.Length > 10)
            {
                double[] logFreq = frequencies.Skip(1).Take(10).Where(f => f > 1e-4).Select(f => Math.Log10(f)).ToArray();
                double[] logPower = psd.Skip(1).Take(10).Where((p, i) => frequencies[i + 1] > 1e-4).Select(p => Math.Log10(Math.Max(p, 1e-10))).ToArray();
                if (logFreq.Length >= 2)
                {
                    double slope = CalculateSlope(logFreq, logPower);
                    if (Math.Abs(slope + 5.0 / 3.0) < 0.5) // Allow some deviation from -5/3
                        turbulenceDetected = true;
                }
            }

            // Output results safely
            if (InvokeRequired)
            {
                Invoke((Action)(() =>
                {
                    if (!outputTextBox.IsDisposed)
                    {
                        outputTextBox.AppendText($"Spectral Analysis Results for {variableComboBox.SelectedItem} (Welch's PSD, {selectedWindowType} Window):\r\n");
                        foreach (var (freq, power, label) in dominantFrequencies)
                        {
                            double period = freq > 0 ? 1.0 / freq : double.MaxValue;
                            outputTextBox.AppendText($"Frequency: {freq:F6} Hz, Period: {period:F2} s, PSD: {power:F4}, Label: {label}\r\n");
                        }

                        outputTextBox.AppendText(turbulenceDetected
                            ? "Turbulence detected: Spectral slope close to -5/3, indicating possible energy cascade.\r\n"
                            : "No clear turbulence signature detected in the spectral slope.\r\n");

                        if (mixingEvents.Count > 0)
                        {
                            outputTextBox.AppendText($"Detected {mixingEvents.Count} potential mixing events:\r\n");
                            foreach (var evt in mixingEvents)
                            {
                                outputTextBox.AppendText($"Time: {evt.Time:F2} s, Value: {evt.Value:F2}\r\n");
                            }
                        }
                        else
                        {
                            outputTextBox.AppendText("No mixing events detected.\r\n");
                        }
                    }
                }));
            }
            else
            {
                if (!outputTextBox.IsDisposed)
                {
                    outputTextBox.AppendText($"Spectral Analysis Results for {variableComboBox.SelectedItem} (Welch's PSD, {selectedWindowType} Window):\r\n");
                    foreach (var (freq, power, label) in dominantFrequencies)
                    {
                        double period = freq > 0 ? 1.0 / freq : double.MaxValue;
                        outputTextBox.AppendText($"Frequency: {freq:F6} Hz, Period: {period:F2} s, PSD: {power:F4}, Label: {label}\r\n");
                    }

                    outputTextBox.AppendText(turbulenceDetected
                        ? "Turbulence detected: Spectral slope close to -5/3, indicating possible energy cascade.\r\n"
                        : "No clear turbulence signature detected in the spectral slope.\r\n");

                    if (mixingEvents.Count > 0)
                    {
                        outputTextBox.AppendText($"Detected {mixingEvents.Count} potential mixing events:\r\n");
                        foreach (var evt in mixingEvents)
                        {
                            outputTextBox.AppendText($"Time: {evt.Time:F2} s, Value: {evt.Value:F2}\r\n");
                        }
                    }
                    else
                    {
                        outputTextBox.AppendText("No mixing events detected.\r\n");
                    }
                }
            }

            if (!visualizationPanel.IsDisposed)
            {
                visualizationPanel.Tag = new SpectralData { PowerSpectrum = psd, Frequencies = frequencies, SegmentLength = segmentLength };
                visualizationPanel.Invalidate();
            }
        }

        private double CalculateSlope(double[] x, double[] y)
        {
            if (x.Length != y.Length || x.Length < 2) return 0;
            double meanX = x.Average();
            double meanY = y.Average();
            double numerator = 0, denominator = 0;
            for (int i = 0; i < x.Length; i++)
            {
                numerator += (x[i] - meanX) * (y[i] - meanY);
                denominator += (x[i] - meanX) * (x[i] - meanX);
            }
            return denominator > 1E-10 ? numerator / denominator : 0;
        }

        private void VisualizationPanel_Paint(object sender, PaintEventArgs e)
        {
            if (IsDisposed || visualizationPanel.IsDisposed) return;

            Graphics g = e.Graphics;
            g.Clear(Color.White);
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

            const int margin = 50;
            float plotWidth = visualizationPanel.Width - 2 * margin;
            float plotHeight = visualizationPanel.Height - 2 * margin - 30; // Extra space for title
            float xOrigin = margin;
            float yOrigin = visualizationPanel.Height - margin;

            bool isLogLog = scaleComboBox.SelectedItem?.ToString() == "Log-Log";

            if (visualizationPanel.Tag is SpectralData spectralData)
            {
                double[] powerSpectrum = spectralData.PowerSpectrum;
                double[] frequencies = spectralData.Frequencies;
                int segmentLength = spectralData.SegmentLength;
                double maxPower = powerSpectrum.Length > 0 ? powerSpectrum.Max() : 1.0;
                double maxFreq = frequencies.Length > 0 ? frequencies[frequencies.Length - 1] : 1.0;
                double minFreq = frequencies.Length > 0 ? frequencies[1] : 1e-6; // Skip DC component

                // Ensure valid values for logarithmic calculations
                if (maxPower < 1e-10) maxPower = 1.0;
                if (minFreq < 1e-6) minFreq = 1e-6;
                if (maxFreq <= minFreq) maxFreq = minFreq + 1e-6;

                // Draw title with window type
                g.DrawString($"Power Spectral Density of {variableComboBox.SelectedItem} ({(isLogLog ? "Log-Log" : "Linear")}, {selectedWindowType} Window)", new Font("Consolas", 10F, FontStyle.Bold), Brushes.Black, margin, 10);

                // Draw axes
                g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin + plotWidth, yOrigin); // X-axis
                g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin, yOrigin - plotHeight); // Y-axis

                // Draw axis labels
                g.DrawString("Frequency (Hz)", new Font("Consolas", 9F), Brushes.Black, xOrigin + plotWidth / 2 - 30, yOrigin + 20);
                g.DrawString("PSD", new Font("Consolas", 9F), Brushes.Black, xOrigin - 40, yOrigin - plotHeight / 2 - 10, new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical });

                // Draw grid and ticks
                DrawGridAndTicks(g, xOrigin, yOrigin, plotWidth, plotHeight, maxFreq, maxPower, true, margin, isLogLog, minFreq);

                // Plot power spectrum
                for (int i = 1; i < powerSpectrum.Length; i++)
                {
                    // Skip invalid frequency or power values
                    if (frequencies[i - 1] <= 0 || frequencies[i] <= 0 || double.IsNaN(powerSpectrum[i - 1]) || double.IsNaN(powerSpectrum[i]))
                        continue;

                    double xScale1 = isLogLog ? (Math.Log10(frequencies[i - 1]) - Math.Log10(minFreq)) / (Math.Log10(maxFreq) - Math.Log10(minFreq)) : frequencies[i - 1] / maxFreq;
                    double xScale2 = isLogLog ? (Math.Log10(frequencies[i]) - Math.Log10(minFreq)) / (Math.Log10(maxFreq) - Math.Log10(minFreq)) : frequencies[i] / maxFreq;
                    double yScale1 = isLogLog ? (Math.Log10(Math.Max(powerSpectrum[i - 1], 1e-10)) - Math.Log10(1e-10)) / (Math.Log10(maxPower) - Math.Log10(1e-10)) : powerSpectrum[i - 1] / maxPower;
                    double yScale2 = isLogLog ? (Math.Log10(Math.Max(powerSpectrum[i], 1e-10)) - Math.Log10(1e-10)) / (Math.Log10(maxPower) - Math.Log10(1e-10)) : powerSpectrum[i] / maxPower;

                    // Check for NaN or Infinity in scaling
                    if (double.IsNaN(xScale1) || double.IsInfinity(xScale1) || double.IsNaN(xScale2) || double.IsInfinity(xScale2) ||
                        double.IsNaN(yScale1) || double.IsInfinity(yScale1) || double.IsNaN(yScale2) || double.IsInfinity(yScale2))
                        continue;

                    float x1 = xOrigin + (float)(xScale1 * plotWidth);
                    float y1 = yOrigin - (float)(yScale1 * plotHeight);
                    float x2 = xOrigin + (float)(xScale2 * plotWidth);
                    float y2 = yOrigin - (float)(yScale2 * plotHeight);

                    // Clamp all coordinates to prevent overflow
                    x1 = Math.Max(margin, Math.Min(xOrigin + plotWidth, x1));
                    y1 = Math.Max(margin, Math.Min(yOrigin, y1));
                    x2 = Math.Max(margin, Math.Min(xOrigin + plotWidth, x2));
                    y2 = Math.Max(margin, Math.Min(yOrigin, y2));

                    g.DrawLine(Pens.Blue, x1, y1, x2, y2);
                }

                // Highlight tidal frequencies (M2, K1)
                double m2Freq = 1.0 / 44712; // M2: ~12.42 hours
                double k1Freq = 1.0 / 86148; // K1: ~23.93 hours
                double freqTolerance = (1.0 / (samplingRate * segmentLength)) * 2;
                foreach (var freq in new[] { (m2Freq, "M2"), (k1Freq, "K1") })
                {
                    if (freq.Item1 <= maxFreq && freq.Item1 >= minFreq)
                    {
                        double xScale = isLogLog ? (Math.Log10(freq.Item1) - Math.Log10(minFreq)) / (Math.Log10(maxFreq) - Math.Log10(minFreq)) : freq.Item1 / maxFreq;
                        if (double.IsNaN(xScale) || double.IsInfinity(xScale))
                            continue;

                        float xPos = xOrigin + (float)(xScale * plotWidth);
                        xPos = Math.Max(margin, Math.Min(xOrigin + plotWidth, xPos));
                        g.DrawLine(Pens.Red, xPos, yOrigin, xPos, yOrigin - plotHeight);
                        g.DrawString(freq.Item2, new Font("Consolas", 8F), Brushes.Red, xPos + 2, yOrigin - plotHeight);
                    }
                }

                // Draw legend
                g.DrawLine(Pens.Blue, xOrigin + plotWidth - 60, 30, xOrigin + plotWidth - 40, 30);
                g.DrawString("Power Spectral Density", new Font("Consolas", 8F), Brushes.Black, xOrigin + plotWidth - 100, 25);
            }
            else if (isAnalyzing)
            {
                double maxValue = timeSeriesData.Count > 0 ? Math.Max(Math.Abs(timeSeriesData.Min()), timeSeriesData.Max()) : 1.0;
                double maxTime = duration;

                // Ensure maxValue is not too small to avoid division issues
                if (maxValue < 1E-10) maxValue = 1.0;

                // Draw title
                g.DrawString($"Time Series: {variableComboBox.SelectedItem?.ToString() ?? "Value"}", new Font("Consolas", 10F, FontStyle.Bold), Brushes.Black, margin, 10);

                // Draw axes
                g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin + plotWidth, yOrigin); // X-axis
                g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin, yOrigin - plotHeight); // Y-axis

                // Draw axis labels
                g.DrawString("Time (s)", new Font("Consolas", 9F), Brushes.Black, xOrigin + plotWidth / 2 - 20, yOrigin + 20);
                g.DrawString(variableComboBox.SelectedItem?.ToString() ?? "Value", new Font("Consolas", 9F), Brushes.Black, xOrigin - 40, yOrigin - plotHeight / 2 - 10, new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical });

                // Draw grid and ticks
                DrawGridAndTicks(g, xOrigin, yOrigin, plotWidth, plotHeight, maxTime, maxValue, false, margin, false, 0);

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

                // Highlight mixing events
                foreach (var evt in mixingEvents)
                {
                    float x = xOrigin + (float)(evt.Time / maxTime * plotWidth);
                    float y = yOrigin - (float)(evt.Value / maxValue * plotHeight / 2 + plotHeight / 2);
                    y = Math.Max(margin, Math.Min(yOrigin, y));
                    g.FillEllipse(Brushes.Red, x - 3, y - 3, 6, 6);
                }

                // Draw legend
                g.DrawLine(Pens.Blue, xOrigin + plotWidth - 60, 30, xOrigin + plotWidth - 40, 30);
                g.DrawString("Time Series", new Font("Consolas", 8F), Brushes.Black, xOrigin + plotWidth - 100, 25);
                if (mixingEvents.Count > 0)
                {
                    g.FillEllipse(Brushes.Red, xOrigin + plotWidth - 60, 40, 6, 6);
                    g.DrawString("Mixing Events", new Font("Consolas", 8F), Brushes.Black, xOrigin + plotWidth - 100, 40);
                }
            }
        }

        private void DrawGridAndTicks(Graphics g, float xOrigin, float yOrigin, float plotWidth, float plotHeight, double maxX, double maxY, bool isPowerSpectrum, int margin, bool isLogLog, double minFreq)
        {
            // Ensure maxX and maxY are not too small to avoid division issues
            if (maxX < 1E-10) maxX = 1.0;
            if (maxY < 1E-10) maxY = 1.0;
            if (isPowerSpectrum && isLogLog && minFreq < 1e-6) minFreq = 1e-6;

            // Calculate tick intervals
            double xTickInterval = isPowerSpectrum && isLogLog ? CalculateLogTickInterval(minFreq, maxX, 5) : CalculateTickInterval(maxX, 5);
            double yTickInterval = isPowerSpectrum && isLogLog ? CalculateLogTickInterval(1e-10, maxY, 5) : CalculateTickInterval(maxY, 5);

            // Draw vertical grid lines and x-axis ticks
            if (isPowerSpectrum && isLogLog)
            {
                double logMinFreq = Math.Log10(minFreq);
                double logMaxFreq = Math.Log10(maxX);
                for (double logX = Math.Ceiling(logMinFreq); logX <= logMaxFreq; logX += xTickInterval)
                {
                    double x = Math.Pow(10, logX);
                    float xPos = xOrigin + (float)((logX - logMinFreq) / (logMaxFreq - logMinFreq) * plotWidth);
                    xPos = Math.Max(margin, Math.Min(xOrigin + plotWidth, xPos));
                    g.DrawLine(Pens.LightGray, xPos, yOrigin, xPos, yOrigin - plotHeight); // Grid line
                    g.DrawLine(Pens.Black, xPos, yOrigin, xPos, yOrigin + 5); // Tick
                    string label = $"10^{logX:F0}";
                    g.DrawString(label, new Font("Consolas", 8F), Brushes.Black, xPos - 15, yOrigin + 10);
                }
            }
            else
            {
                for (double x = 0; x <= maxX; x += xTickInterval)
                {
                    float xPos = xOrigin + (float)(x / maxX * plotWidth);
                    xPos = Math.Max(margin, Math.Min(xOrigin + plotWidth, xPos));
                    g.DrawLine(Pens.LightGray, xPos, yOrigin, xPos, yOrigin - plotHeight); // Grid line
                    g.DrawLine(Pens.Black, xPos, yOrigin, xPos, yOrigin + 5); // Tick
                    string label = isPowerSpectrum ? x.ToString("F4") : ((int)x).ToString();
                    g.DrawString(label, new Font("Consolas", 8F), Brushes.Black, xPos - 15, yOrigin + 10);
                }
            }

            // Draw horizontal grid lines and y-axis ticks
            double yStart = isPowerSpectrum ? 0 : -maxY;
            if (isPowerSpectrum && isLogLog)
            {
                double logMinPower = Math.Log10(1e-10);
                double logMaxPower = Math.Log10(maxY);
                for (double logY = Math.Ceiling(logMinPower); logY <= logMaxPower; logY += yTickInterval)
                {
                    double y = Math.Pow(10, logY);
                    float yPos = yOrigin - (float)((logY - logMinPower) / (logMaxPower - logMinPower) * plotHeight);
                    yPos = Math.Max(margin, Math.Min(yOrigin, yPos));
                    g.DrawLine(Pens.LightGray, xOrigin, yPos, xOrigin + plotWidth, yPos); // Grid line
                    g.DrawLine(Pens.Black, xOrigin - 5, yPos, xOrigin, yPos); // Tick
                    string label = $"10^{logY:F0}";
                    g.DrawString(label, new Font("Consolas", 8F), Brushes.Black, xOrigin - 40, yPos - 5);
                }
            }
            else
            {
                for (double y = yStart; y <= maxY; y += yTickInterval)
                {
                    float yPos = yOrigin - (float)(y / maxY * plotHeight / (isPowerSpectrum ? 1 : 2) + (isPowerSpectrum ? 0 : plotHeight / 2));
                    yPos = Math.Max(margin, Math.Min(yOrigin, yPos));
                    g.DrawLine(Pens.LightGray, xOrigin, yPos, xOrigin + plotWidth, yPos); // Grid line
                    g.DrawLine(Pens.Black, xOrigin - 5, yPos, xOrigin, yPos); // Tick
                    string label = y.ToString("F2");
                    g.DrawString(label, new Font("Consolas", 8F), Brushes.Black, xOrigin - 40, yPos - 5);
                }
            }
        }

        private double CalculateTickInterval(double range, int targetTicks)
        {
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

        private double CalculateLogTickInterval(double minVal, double maxVal, int targetTicks)
        {
            double logMin = Math.Log10(minVal);
            double logMax = Math.Log10(maxVal);
            return (logMax - logMin) / targetTicks;
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