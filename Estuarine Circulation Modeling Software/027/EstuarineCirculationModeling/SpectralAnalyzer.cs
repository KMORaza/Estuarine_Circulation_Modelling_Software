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
        private ComboBox analysisComboBox;
        private ComboBox modeComboBox;
        private TextBox durationTextBox;
        private TextBox samplingRateTextBox;
        private Button analyzeButton;
        private Button pauseButton;
        private Button resetButton;
        private Button clearButton;
        private Panel visualizationPanel;
        private Panel heatmapPanel;
        private TextBox outputTextBox;
        private Timer simulationTimer;
        private List<double> timeSeriesData;
        private List<double[,]> salinityData; // [space, time]
        private List<double[,]> velocityData; // [depth, time]
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

        // Custom class to hold heatmap data
        private class HeatmapData
        {
            public double[,] PowerSpectrum { get; set; } // [frequency, time]
            public double[] Frequencies { get; set; }
            public double[] Times { get; set; }
            public int SegmentLength { get; set; }
        }

        // Custom class to hold EOF/POD data
        private class EOFData
        {
            public double[][] SpatialModes { get; set; } // [mode][space]
            public double[][] TemporalCoefficients { get; set; } // [mode][time]
            public double[] SingularValues { get; set; }
            public int NumModes { get; set; }
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
            salinityData = new List<double[,]>();
            velocityData = new List<double[,]>();
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
            this.Size = new Size(1200, 600);
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

            // Analysis Selection
            Label analysisLabel = new Label
            {
                Text = "Analysis Type:",
                Location = new Point(10, 10),
                AutoSize = true
            };

            analysisComboBox = new ComboBox
            {
                Location = new Point(10, 30),
                Size = new Size(200, 22),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            analysisComboBox.Items.AddRange(new string[] { "Spectral Analysis", "EOF/POD Analysis" });
            analysisComboBox.SelectedIndex = 0;
            analysisComboBox.SelectedIndexChanged += (s, e) =>
            {
                if (!IsDisposed)
                {
                    bool isSpectral = analysisComboBox.SelectedItem?.ToString() == "Spectral Analysis";
                    variableComboBox.Items.Clear();
                    variableComboBox.Items.AddRange(isSpectral ?
                        new string[] { "Richardson Number", "Velocity", "Water Level" } :
                        new string[] { "Salinity", "Velocity Profile" });
                    variableComboBox.SelectedIndex = 0;
                    windowComboBox.Enabled = isSpectral;
                    modeComboBox.Visible = !isSpectral;
                    visualizationPanel.Invalidate();
                    heatmapPanel.Invalidate();
                }
            };

            // Variable Selection
            Label variableLabel = new Label
            {
                Text = "Select Variable:",
                Location = new Point(10, 60),
                AutoSize = true
            };

            variableComboBox = new ComboBox
            {
                Location = new Point(10, 80),
                Size = new Size(200, 22),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            variableComboBox.Items.AddRange(new string[] { "Richardson Number", "Velocity", "Water Level" });
            variableComboBox.SelectedIndex = 0;

            // Scale Selection
            Label scaleLabel = new Label
            {
                Text = "Plot Scale:",
                Location = new Point(10, 110),
                AutoSize = true
            };

            scaleComboBox = new ComboBox
            {
                Location = new Point(10, 130),
                Size = new Size(200, 22),
                DropDownStyle = ComboBoxStyle.DropDownList
            };
            scaleComboBox.Items.AddRange(new string[] { "Linear", "Log-Log" });
            scaleComboBox.SelectedIndex = 0;
            scaleComboBox.SelectedIndexChanged += (s, e) =>
            {
                if (!IsDisposed)
                {
                    visualizationPanel.Invalidate();
                    heatmapPanel.Invalidate();
                }
            };

            // Window Selection
            Label windowLabel = new Label
            {
                Text = "Window Type:",
                Location = new Point(10, 160),
                AutoSize = true
            };

            windowComboBox = new ComboBox
            {
                Location = new Point(10, 180),
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
                    visualizationPanel.Invalidate();
                    heatmapPanel.Invalidate();
                }
            };

            // Mode Selection for EOF/POD
            Label modeLabel = new Label
            {
                Text = "EOF/POD Mode:",
                Location = new Point(10, 210),
                AutoSize = true
            };

            modeComboBox = new ComboBox
            {
                Location = new Point(10, 230),
                Size = new Size(200, 22),
                DropDownStyle = ComboBoxStyle.DropDownList,
                Visible = false
            };
            modeComboBox.Items.AddRange(new string[] { "Mode 1", "Mode 2", "Mode 3" });
            modeComboBox.SelectedIndex = 0;
            modeComboBox.SelectedIndexChanged += (s, e) =>
            {
                if (!IsDisposed)
                {
                    visualizationPanel.Invalidate();
                    heatmapPanel.Invalidate();
                }
            };

            // Duration
            Label durationLabel = new Label
            {
                Text = "Duration (s):",
                Location = new Point(10, 260),
                AutoSize = true
            };

            durationTextBox = new TextBox
            {
                Location = new Point(10, 280),
                Size = new Size(200, 22),
                Text = "86400"
            };

            // Sampling Rate
            Label samplingRateLabel = new Label
            {
                Text = "Sampling Rate (s):",
                Location = new Point(10, 310),
                AutoSize = true
            };

            samplingRateTextBox = new TextBox
            {
                Location = new Point(10, 330),
                Size = new Size(200, 22),
                Text = "3600"
            };

            // Buttons
            analyzeButton = new Button
            {
                Text = "Analyze",
                Location = new Point(10, 360),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            analyzeButton.FlatAppearance.BorderSize = 1;
            analyzeButton.FlatAppearance.BorderColor = Color.Black;
            analyzeButton.Click += AnalyzeButton_Click;

            pauseButton = new Button
            {
                Text = "Pause",
                Location = new Point(120, 360),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            pauseButton.FlatAppearance.BorderSize = 1;
            pauseButton.FlatAppearance.BorderColor = Color.Black;
            pauseButton.Click += PauseButton_Click;

            resetButton = new Button
            {
                Text = "Reset",
                Location = new Point(10, 395),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            resetButton.FlatAppearance.BorderSize = 1;
            resetButton.FlatAppearance.BorderColor = Color.Black;
            resetButton.Click += ResetButton_Click;

            clearButton = new Button
            {
                Text = "Clear",
                Location = new Point(120, 395),
                Size = new Size(100, 25),
                FlatStyle = FlatStyle.Flat
            };
            clearButton.FlatAppearance.BorderSize = 1;
            clearButton.FlatAppearance.BorderColor = Color.Black;
            clearButton.Click += ClearButton_Click;

            // Visualization Panel (PSD/Time Series or EOF Spatial Modes)
            visualizationPanel = new Panel
            {
                Location = new Point(270, 10),
                Size = new Size(450, 350),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            visualizationPanel.Paint += VisualizationPanel_Paint;

            // Heatmap Panel (Spectral Heatmap or EOF Temporal Coefficients)
            heatmapPanel = new Panel
            {
                Location = new Point(730, 10),
                Size = new Size(450, 350),
                BorderStyle = BorderStyle.FixedSingle,
                BackColor = Color.White
            };
            heatmapPanel.Paint += HeatmapPanel_Paint;

            // Output TextBox
            outputTextBox = new TextBox
            {
                Location = new Point(270, 370),
                Size = new Size(910, 150),
                Multiline = true,
                ReadOnly = true,
                ScrollBars = ScrollBars.Vertical,
                Font = new Font("Consolas", 9F)
            };

            // Add controls to control panel
            controlPanel.Controls.AddRange(new Control[] { analysisLabel, analysisComboBox, variableLabel, variableComboBox, scaleLabel, scaleComboBox, windowLabel, windowComboBox, modeLabel, modeComboBox, durationLabel, durationTextBox, samplingRateLabel, samplingRateTextBox, analyzeButton, pauseButton, resetButton, clearButton });

            // Add controls to form
            this.Controls.AddRange(new Control[] { controlPanel, visualizationPanel, heatmapPanel, outputTextBox });
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
            resetButton.Enabled = timeSeriesData.Count > 0 || salinityData.Count > 0 || velocityData.Count > 0;
            clearButton.Enabled = timeSeriesData.Count > 0 || salinityData.Count > 0 || velocityData.Count > 0 || visualizationPanel.Tag != null || heatmapPanel.Tag != null;
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
                    salinityData.Clear();
                    velocityData.Clear();
                    mixingEvents.Clear();
                    currentTime = 0;
                }
                isAnalyzing = true;
                simulationTimer.Start();
                if (InvokeRequired)
                {
                    Invoke((Action)(() => outputTextBox.AppendText(isAnalyzing && currentTime == 0 ? "Starting analysis simulation...\r\n" : "Resuming analysis simulation...\r\n")));
                }
                else
                {
                    outputTextBox.AppendText(isAnalyzing && currentTime == 0 ? "Starting analysis simulation...\r\n" : "Resuming analysis simulation...\r\n");
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
                Invoke((Action)(() => outputTextBox.AppendText("Analysis paused.\r\n")));
            }
            else
            {
                outputTextBox.AppendText("Analysis paused.\r\n");
            }
            UpdateButtonStates();
        }

        private void ResetButton_Click(object sender, EventArgs e)
        {
            if (IsDisposed) return;

            simulationTimer.Stop();
            isAnalyzing = false;
            timeSeriesData.Clear();
            salinityData.Clear();
            velocityData.Clear();
            mixingEvents.Clear();
            currentTime = 0;
            visualizationPanel.Tag = null;
            heatmapPanel.Tag = null;
            if (!visualizationPanel.IsDisposed) visualizationPanel.Invalidate();
            if (!heatmapPanel.IsDisposed) heatmapPanel.Invalidate();
            if (InvokeRequired)
            {
                Invoke((Action)(() => outputTextBox.AppendText("Analysis reset.\r\n")));
            }
            else
            {
                outputTextBox.AppendText("Analysis reset.\r\n");
            }
            UpdateButtonStates();
        }

        private void ClearButton_Click(object sender, EventArgs e)
        {
            if (IsDisposed) return;

            simulationTimer.Stop();
            isAnalyzing = false;
            timeSeriesData.Clear();
            salinityData.Clear();
            velocityData.Clear();
            mixingEvents.Clear();
            currentTime = 0;
            if (InvokeRequired)
            {
                Invoke((Action)(() => outputTextBox.Clear()));
                Invoke((Action)(() => outputTextBox.AppendText("Analysis cleared.\r\n")));
            }
            else
            {
                outputTextBox.Clear();
                outputTextBox.AppendText("Analysis cleared.\r\n");
            }
            visualizationPanel.Tag = null;
            heatmapPanel.Tag = null;
            if (!visualizationPanel.IsDisposed) visualizationPanel.Invalidate();
            if (!heatmapPanel.IsDisposed) heatmapPanel.Invalidate();
            UpdateButtonStates();
        }

        private void UpdateSimulation()
        {
            if (IsDisposed) return;

            if (currentTime >= duration)
            {
                simulationTimer.Stop();
                isAnalyzing = false;
                PerformAnalysis();
                UpdateButtonStates();
                return;
            }

            // Simulate data based on analysis type
            string analysisType = analysisComboBox.SelectedItem?.ToString() ?? "Spectral Analysis";
            string selectedVariable = variableComboBox.SelectedItem?.ToString() ?? "Unknown";

            if (analysisType == "Spectral Analysis")
            {
                double value = GenerateSyntheticData(selectedVariable, currentTime);
                timeSeriesData.Add(value);

                // Detect mixing events (high variability)
                if (timeSeriesData.Count >= 5)
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
            }
            else // EOF/POD Analysis
            {
                if (selectedVariable == "Salinity")
                {
                    double[,] salinity = GenerateSyntheticSalinityData(currentTime);
                    salinityData.Add(salinity);
                }
                else // Velocity Profile
                {
                    double[,] velocity = GenerateSyntheticVelocityData(currentTime);
                    velocityData.Add(velocity);
                }
            }

            currentTime += samplingRate;
            if (!visualizationPanel.IsDisposed) visualizationPanel.Invalidate();
            UpdateButtonStates();
        }

        private double GenerateSyntheticData(string variable, double time)
        {
            double tidalPeriod = 43200; // 12-hour tidal cycle
            double inertialPeriod = 17 * 3600; // Approx inertial period at 45° latitude
            double value = 0;

            if (variable == "Richardson Number")
            {
                value = 0.5 + 0.3 * Math.Sin(2 * Math.PI * time / tidalPeriod) + 0.2 * Math.Sin(2 * Math.PI * time / inertialPeriod);
                if (new Random((int)time).NextDouble() < 0.02) value += 0.5;
            }
            else if (variable == "Velocity")
            {
                value = 0.1 + 0.05 * Math.Sin(2 * Math.PI * time / tidalPeriod) + 0.02 * Math.Sin(2 * Math.PI * time / (tidalPeriod / 2)) + new Random((int)time).NextDouble() * 0.01;
            }
            else if (variable == "Water Level")
            {
                value = 1.0 * Math.Sin(2 * Math.PI * time / 44712) + 0.5 * Math.Sin(2 * Math.PI * time / 86148); // M2 + K1
            }

            return value;
        }

        private double[,] GenerateSyntheticSalinityData(double time)
        {
            const int gridSizeX = 10;
            const int gridSizeY = 10;
            double[,] salinity = new double[gridSizeX * gridSizeY, 1];
            double tidalPeriod = 44712; // M2 tide
            Random rand = new Random((int)time);
            int idx = 0;
            for (int x = 0; x < gridSizeX; x++)
            {
                for (int y = 0; y < gridSizeY; y++)
                {
                    double spatialFactor = 1.0 - 0.1 * (x + y); // Gradient across grid
                    double tidal = 30.0 + 2.0 * Math.Sin(2 * Math.PI * time / tidalPeriod);
                    double noise = rand.NextDouble() * 0.5;
                    salinity[idx++, 0] = tidal * spatialFactor + noise;
                }
            }
            return salinity;
        }

        private double[,] GenerateSyntheticVelocityData(double time)
        {
            const int depthLevels = 10;
            double[,] velocity = new double[depthLevels, 1];
            double tidalPeriod = 44712; // M2 tide
            Random rand = new Random((int)time);
            for (int z = 0; z < depthLevels; z++)
            {
                double depthFactor = 1.0 - 0.1 * z; // Decrease with depth
                double tidal = 0.1 * Math.Sin(2 * Math.PI * time / tidalPeriod);
                double noise = rand.NextDouble() * 0.01;
                velocity[z, 0] = tidal * depthFactor + noise;
            }
            return velocity;
        }

        private void PerformAnalysis()
        {
            if (IsDisposed) return;

            string analysisType = analysisComboBox.SelectedItem?.ToString() ?? "Spectral Analysis";
            if (analysisType == "Spectral Analysis")
            {
                PerformSpectralAnalysis();
            }
            else
            {
                PerformEOFAnalysis();
            }
        }

        private void PerformSpectralAnalysis()
        {
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
            int maxSegmentLength = Math.Min(256, data.Length);
            // Ensure segmentLength is a power of 2
            int segmentLength = (int)Math.Pow(2, Math.Floor(Math.Log(maxSegmentLength, 2)));
            if (segmentLength < 2) segmentLength = 2; // Minimum power of 2
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

            // Compute STFT for heatmap
            double[,] heatmapPsd = new double[fftLength, numSegments];
            double[] heatmapTimes = new double[numSegments];
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
                    double power = segmentData[i].Magnitude * segmentData[i].Magnitude / (windowPower * samplingRate);
                    psd[i] += power; // For Welch's PSD
                    heatmapPsd[i, seg] = power; // For heatmap
                }
                heatmapTimes[seg] = start * samplingRate;
            }

            // Average PSD for Welch's method
            for (int i = 0; i < fftLength; i++)
            {
                psd[i] /= numSegments;
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

            // Analyze dominant frequency evolution for heatmap
            List<(double Time, double Frequency, double Power)> dominantFreqEvolution = new List<(double, double, double)>();
            for (int seg = 0; seg < numSegments; seg++)
            {
                double maxPower = 0;
                int maxIndex = 1; // Skip DC
                for (int i = 1; i < fftLength; i++)
                {
                    if (heatmapPsd[i, seg] > maxPower)
                    {
                        maxPower = heatmapPsd[i, seg];
                        maxIndex = i;
                    }
                }
                dominantFreqEvolution.Add((heatmapTimes[seg], frequencies[maxIndex], maxPower));
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

                        outputTextBox.AppendText($"Time-Frequency Heatmap Analysis ({selectedWindowType} Window):\r\n");
                        foreach (var (time, freq, power) in dominantFreqEvolution)
                        {
                            double period = freq > 0 ? 1.0 / freq : double.MaxValue;
                            outputTextBox.AppendText($"Time: {time:F2} s, Dominant Frequency: {freq:F6} Hz, Period: {period:F2} s, Power: {power:F4}\r\n");
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

                    outputTextBox.AppendText($"Time-Frequency Heatmap Analysis ({selectedWindowType} Window):\r\n");
                    foreach (var (time, freq, power) in dominantFreqEvolution)
                    {
                        double period = freq > 0 ? 1.0 / freq : double.MaxValue;
                        outputTextBox.AppendText($"Time: {time:F2} s, Dominant Frequency: {freq:F6} Hz, Period: {period:F2} s, Power: {power:F4}\r\n");
                    }
                }
            }

            if (!visualizationPanel.IsDisposed)
            {
                visualizationPanel.Tag = new SpectralData { PowerSpectrum = psd, Frequencies = frequencies, SegmentLength = segmentLength };
                visualizationPanel.Invalidate();
            }
            if (!heatmapPanel.IsDisposed)
            {
                heatmapPanel.Tag = new HeatmapData
                {
                    PowerSpectrum = heatmapPsd,
                    Frequencies = frequencies,
                    Times = heatmapTimes,
                    SegmentLength = segmentLength
                };
                heatmapPanel.Invalidate();
            }
        }

        private void PerformEOFAnalysis()
        {
            if (IsDisposed) return;

            string selectedVariable = variableComboBox.SelectedItem?.ToString() ?? "Salinity";
            double[,] dataMatrix;
            int spatialDim, timeDim;

            if (selectedVariable == "Salinity")
            {
                if (salinityData.Count < 2)
                {
                    if (InvokeRequired)
                        Invoke((Action)(() => outputTextBox.AppendText("Insufficient salinity data for EOF analysis.\r\n")));
                    else
                        outputTextBox.AppendText("Insufficient salinity data for EOF analysis.\r\n");
                    return;
                }
                spatialDim = salinityData[0].GetLength(0); // 100 (10x10 grid)
                timeDim = salinityData.Count;
                dataMatrix = new double[spatialDim, timeDim];
                for (int t = 0; t < timeDim; t++)
                    for (int s = 0; s < spatialDim; s++)
                        dataMatrix[s, t] = salinityData[t][s, 0];
            }
            else // Velocity Profile
            {
                if (velocityData.Count < 2)
                {
                    if (InvokeRequired)
                        Invoke((Action)(() => outputTextBox.AppendText("Insufficient velocity data for EOF analysis.\r\n")));
                    else
                        outputTextBox.AppendText("Insufficient velocity data for EOF analysis.\r\n");
                    return;
                }
                spatialDim = velocityData[0].GetLength(0); // 10 (depth levels)
                timeDim = velocityData.Count;
                dataMatrix = new double[spatialDim, timeDim];
                for (int t = 0; t < timeDim; t++)
                    for (int s = 0; s < spatialDim; s++)
                        dataMatrix[s, t] = velocityData[t][s, 0];
            }

            // Subtract mean from each spatial point
            for (int s = 0; s < spatialDim; s++)
            {
                double mean = 0;
                for (int t = 0; t < timeDim; t++)
                    mean += dataMatrix[s, t];
                mean /= timeDim;
                for (int t = 0; t < timeDim; t++)
                    dataMatrix[s, t] -= mean;
            }

            // Perform SVD (simplified power iteration for top 3 modes)
            int numModes = Math.Min(3, Math.Min(spatialDim, timeDim));
            double[][] spatialModes = new double[numModes][];
            double[][] temporalCoefficients = new double[numModes][];
            double[] singularValues = new double[numModes];

            for (int m = 0; m < numModes; m++)
            {
                // Initialize random vector
                double[] v = new double[timeDim];
                Random rand = new Random(m);
                for (int i = 0; i < timeDim; i++)
                    v[i] = rand.NextDouble();
                NormalizeVector(v);

                // Power iteration
                for (int iter = 0; iter < 20; iter++)
                {
                    double[] u = new double[spatialDim];
                    // Compute A * v
                    for (int i = 0; i < spatialDim; i++)
                    {
                        u[i] = 0;
                        for (int j = 0; j < timeDim; j++)
                            u[i] += dataMatrix[i, j] * v[j];
                    }
                    NormalizeVector(u);

                    // Compute A^T * u
                    double[] vNew = new double[timeDim];
                    for (int j = 0; j < timeDim; j++)
                    {
                        vNew[j] = 0;
                        for (int i = 0; i < spatialDim; i++)
                            vNew[j] += dataMatrix[i, j] * u[i];
                    }
                    NormalizeVector(vNew);
                    v = vNew;
                }

                // Compute singular value
                double[] Av = new double[spatialDim];
                for (int i = 0; i < spatialDim; i++)
                {
                    Av[i] = 0;
                    for (int j = 0; j < timeDim; j++)
                        Av[i] += dataMatrix[i, j] * v[j];
                }
                double sigma = Math.Sqrt(Av.Sum(x => x * x));
                if (sigma < 1e-10) sigma = 1e-10;

                spatialModes[m] = Av.Select(x => x / sigma).ToArray();
                temporalCoefficients[m] = v;
                singularValues[m] = sigma;

                // Deflate matrix
                for (int i = 0; i < spatialDim; i++)
                    for (int j = 0; j < timeDim; j++)
                        dataMatrix[i, j] -= sigma * spatialModes[m][i] * temporalCoefficients[m][j];
            }

            // Calculate explained variance
            double totalVariance = singularValues.Sum(s => s * s);
            double[] explainedVariance = singularValues.Select(s => s * s / totalVariance * 100).ToArray();

            // Output results
            if (InvokeRequired)
            {
                Invoke((Action)(() =>
                {
                    if (!outputTextBox.IsDisposed)
                    {
                        outputTextBox.AppendText($"EOF/POD Analysis Results for {selectedVariable}:\r\n");
                        for (int m = 0; m < numModes; m++)
                        {
                            outputTextBox.AppendText($"Mode {m + 1}: Explained Variance = {explainedVariance[m]:F2}%\r\n");
                        }
                    }
                }));
            }
            else
            {
                if (!outputTextBox.IsDisposed)
                {
                    outputTextBox.AppendText($"EOF/POD Analysis Results for {selectedVariable}:\r\n");
                    for (int m = 0; m < numModes; m++)
                    {
                        outputTextBox.AppendText($"Mode {m + 1}: Explained Variance = {explainedVariance[m]:F2}%\r\n");
                    }
                }
            }

            if (!visualizationPanel.IsDisposed && !heatmapPanel.IsDisposed)
            {
                visualizationPanel.Tag = new EOFData
                {
                    SpatialModes = spatialModes,
                    TemporalCoefficients = temporalCoefficients,
                    SingularValues = singularValues,
                    NumModes = numModes
                };
                heatmapPanel.Tag = visualizationPanel.Tag;
                visualizationPanel.Invalidate();
                heatmapPanel.Invalidate();
            }
        }

        private void NormalizeVector(double[] v)
        {
            double norm = Math.Sqrt(v.Sum(x => x * x));
            if (norm < 1e-10) norm = 1e-10;
            for (int i = 0; i < v.Length; i++)
                v[i] /= norm;
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
            float plotHeight = visualizationPanel.Height - 2 * margin - 30;
            float xOrigin = margin;
            float yOrigin = visualizationPanel.Height - margin;

            bool isLogLog = scaleComboBox.SelectedItem?.ToString() == "Log-Log";
            string analysisType = analysisComboBox.SelectedItem?.ToString() ?? "Spectral Analysis";

            if (analysisType == "Spectral Analysis")
            {
                if (visualizationPanel.Tag is SpectralData spectralData)
                {
                    double[] powerSpectrum = spectralData.PowerSpectrum;
                    double[] frequencies = spectralData.Frequencies;
                    int segmentLength = spectralData.SegmentLength;
                    double maxPower = powerSpectrum.Length > 0 ? powerSpectrum.Max() : 1.0;
                    double maxFreq = frequencies.Length > 0 ? frequencies[frequencies.Length - 1] : 1.0;
                    double minFreq = frequencies.Length > 0 ? frequencies[1] : 1e-6;

                    if (maxPower < 1e-10) maxPower = 1.0;
                    if (minFreq < 1e-6) minFreq = 1e-6;
                    if (maxFreq <= minFreq) maxFreq = minFreq + 1e-6;

                    g.DrawString($"Power Spectral Density of {variableComboBox.SelectedItem} ({(isLogLog ? "Log-Log" : "Linear")}, {selectedWindowType} Window)", new Font("Consolas", 10F, FontStyle.Bold), Brushes.Black, margin, 10);
                    g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin + plotWidth, yOrigin);
                    g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin, yOrigin - plotHeight);
                    g.DrawString("Frequency (Hz)", new Font("Consolas", 9F), Brushes.Black, xOrigin + plotWidth / 2 - 30, yOrigin + 20);
                    g.DrawString("PSD", new Font("Consolas", 9F), Brushes.Black, xOrigin - 40, yOrigin - plotHeight / 2 - 10, new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical });

                    DrawGridAndTicks(g, xOrigin, yOrigin, plotWidth, plotHeight, maxFreq, maxPower, true, margin, isLogLog, minFreq);

                    for (int i = 1; i < powerSpectrum.Length; i++)
                    {
                        if (frequencies[i - 1] <= 0 || frequencies[i] <= 0 || double.IsNaN(powerSpectrum[i - 1]) || double.IsNaN(powerSpectrum[i]))
                            continue;

                        double xScale1 = isLogLog ? (Math.Log10(frequencies[i - 1]) - Math.Log10(minFreq)) / (Math.Log10(maxFreq) - Math.Log10(minFreq)) : frequencies[i - 1] / maxFreq;
                        double xScale2 = isLogLog ? (Math.Log10(frequencies[i]) - Math.Log10(minFreq)) / (Math.Log10(maxFreq) - Math.Log10(minFreq)) : frequencies[i] / maxFreq;
                        double yScale1 = isLogLog ? (Math.Log10(Math.Max(powerSpectrum[i - 1], 1e-10)) - Math.Log10(1e-10)) / (Math.Log10(maxPower) - Math.Log10(1e-10)) : powerSpectrum[i - 1] / maxPower;
                        double yScale2 = isLogLog ? (Math.Log10(Math.Max(powerSpectrum[i], 1e-10)) - Math.Log10(1e-10)) / (Math.Log10(maxPower) - Math.Log10(1e-10)) : powerSpectrum[i] / maxPower;

                        if (double.IsNaN(xScale1) || double.IsInfinity(xScale1) || double.IsNaN(xScale2) || double.IsInfinity(xScale2) ||
                            double.IsNaN(yScale1) || double.IsInfinity(yScale1) || double.IsNaN(yScale2) || double.IsInfinity(yScale2))
                            continue;

                        float x1 = xOrigin + (float)(xScale1 * plotWidth);
                        float y1 = yOrigin - (float)(yScale1 * plotHeight);
                        float x2 = xOrigin + (float)(xScale2 * plotWidth);
                        float y2 = yOrigin - (float)(yScale2 * plotHeight);

                        x1 = Math.Max(margin, Math.Min(xOrigin + plotWidth, x1));
                        y1 = Math.Max(margin, Math.Min(yOrigin, y1));
                        x2 = Math.Max(margin, Math.Min(xOrigin + plotWidth, x2));
                        y2 = Math.Max(margin, Math.Min(yOrigin, y2));

                        g.DrawLine(Pens.Blue, x1, y1, x2, y2);
                    }

                    double m2Freq = 1.0 / 44712;
                    double k1Freq = 1.0 / 86148;
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

                    g.DrawLine(Pens.Blue, xOrigin + plotWidth - 60, 30, xOrigin + plotWidth - 40, 30);
                    g.DrawString("Power Spectral Density", new Font("Consolas", 8F), Brushes.Black, xOrigin + plotWidth - 100, 25);
                }
                else if (isAnalyzing)
                {
                    double maxValue = timeSeriesData.Count > 0 ? Math.Max(Math.Abs(timeSeriesData.Min()), timeSeriesData.Max()) : 1.0;
                    double maxTime = duration;

                    if (maxValue < 1E-10) maxValue = 1.0;

                    g.DrawString($"Time Series: {variableComboBox.SelectedItem?.ToString() ?? "Value"}", new Font("Consolas", 10F, FontStyle.Bold), Brushes.Black, margin, 10);
                    g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin + plotWidth, yOrigin);
                    g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin, yOrigin - plotHeight);
                    g.DrawString("Time (s)", new Font("Consolas", 9F), Brushes.Black, xOrigin + plotWidth / 2 - 20, yOrigin + 20);
                    g.DrawString(variableComboBox.SelectedItem?.ToString() ?? "Value", new Font("Consolas", 9F), Brushes.Black, xOrigin - 40, yOrigin - plotHeight / 2 - 10, new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical });

                    DrawGridAndTicks(g, xOrigin, yOrigin, plotWidth, plotHeight, maxTime, maxValue, false, margin, false, 0);

                    for (int i = 1; i < timeSeriesData.Count; i++)
                    {
                        float x1 = xOrigin + (float)((i - 1) * samplingRate / maxTime * plotWidth);
                        float y1 = yOrigin - (float)(timeSeriesData[i - 1] / maxValue * plotHeight / 2 + plotHeight / 2);
                        float x2 = xOrigin + (float)(i * samplingRate / maxTime * plotWidth);
                        float y2 = yOrigin - (float)(timeSeriesData[i] / maxValue * plotHeight / 2 + plotHeight / 2);
                        y1 = Math.Max(margin, Math.Min(yOrigin, y1));
                        y2 = Math.Max(margin, Math.Min(yOrigin, y2));
                        g.DrawLine(Pens.Blue, x1, y1, x2, y2);
                    }

                    foreach (var evt in mixingEvents)
                    {
                        float x = xOrigin + (float)(evt.Time / maxTime * plotWidth);
                        float y = yOrigin - (float)(evt.Value / maxValue * plotHeight / 2 + plotHeight / 2);
                        y = Math.Max(margin, Math.Min(yOrigin, y));
                        g.FillEllipse(Brushes.Red, x - 3, y - 3, 6, 6);
                    }

                    g.DrawLine(Pens.Blue, xOrigin + plotWidth - 60, 30, xOrigin + plotWidth - 40, 30);
                    g.DrawString("Time Series", new Font("Consolas", 8F), Brushes.Black, xOrigin + plotWidth - 100, 25);
                    if (mixingEvents.Count > 0)
                    {
                        g.FillEllipse(Brushes.Red, xOrigin + plotWidth - 60, 40, 6, 6);
                        g.DrawString("Mixing Events", new Font("Consolas", 8F), Brushes.Black, xOrigin + plotWidth - 100, 40);
                    }
                }
            }
            else // EOF/POD Analysis
            {
                if (visualizationPanel.Tag is EOFData eofData)
                {
                    int modeIndex = modeComboBox.SelectedIndex;
                    if (modeIndex < 0 || modeIndex >= eofData.NumModes) modeIndex = 0;
                    double[] spatialMode = eofData.SpatialModes[modeIndex];
                    string variable = variableComboBox.SelectedItem?.ToString() ?? "Salinity";

                    if (variable == "Salinity")
                    {
                        const int gridSizeX = 10;
                        const int gridSizeY = 10;
                        double maxValue = spatialMode.Length > 0 ? Math.Max(Math.Abs(spatialMode.Min()), Math.Abs(spatialMode.Max())) : 1.0;
                        if (maxValue < 1e-10) maxValue = 1.0;

                        g.DrawString($"EOF Mode {modeIndex + 1} for Salinity", new Font("Consolas", 10F, FontStyle.Bold), Brushes.Black, margin, 10);
                        g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin + plotWidth, yOrigin);
                        g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin, yOrigin - plotHeight);
                        g.DrawString("X", new Font("Consolas", 9F), Brushes.Black, xOrigin + plotWidth / 2 - 10, yOrigin + 20);
                        g.DrawString("Y", new Font("Consolas", 9F), Brushes.Black, xOrigin - 30, yOrigin - plotHeight / 2 - 10, new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical });

                        float cellWidth = plotWidth / gridSizeX;
                        float cellHeight = plotHeight / gridSizeY;
                        for (int x = 0; x < gridSizeX; x++)
                        {
                            for (int y = 0; y < gridSizeY; y++)
                            {
                                int idx = x * gridSizeY + y;
                                double value = spatialMode[idx];
                                double valueScale = value / maxValue;
                                valueScale = Math.Max(-1, Math.Min(1, valueScale));
                                int r = valueScale > 0 ? (int)(valueScale * 255) : 0;
                                int b = valueScale < 0 ? (int)(-valueScale * 255) : 0;
                                using (SolidBrush brush = new SolidBrush(Color.FromArgb(r, 0, b)))
                                {
                                    float px = xOrigin + x * cellWidth;
                                    float py = yOrigin - (y + 1) * cellHeight;
                                    px = Math.Max(margin, Math.Min(xOrigin + plotWidth, px));
                                    py = Math.Max(margin, Math.Min(yOrigin, py));
                                    g.FillRectangle(brush, px, py, cellWidth, cellHeight);
                                }
                            }
                        }

                        // Draw grid
                        for (int x = 0; x <= gridSizeX; x++)
                        {
                            float px = xOrigin + x * cellWidth;
                            g.DrawLine(Pens.LightGray, px, yOrigin, px, yOrigin - plotHeight);
                        }
                        for (int y = 0; y <= gridSizeY; y++)
                        {
                            float py = yOrigin - y * cellHeight;
                            g.DrawLine(Pens.LightGray, xOrigin, py, xOrigin + plotWidth, py);
                        }

                        // Draw legend
                        float legendWidth = 20;
                        float legendHeight = plotHeight / 2;
                        float legendX = xOrigin + plotWidth - 30;
                        float legendY = yOrigin - legendHeight;
                        for (int i = 0; i < (int)legendHeight; i++)
                        {
                            float valueScale = 1 - 2 * i / legendHeight;
                            int r = valueScale > 0 ? (int)(valueScale * 255) : 0;
                            int b = valueScale < 0 ? (int)(-valueScale * 255) : 0;
                            using (SolidBrush brush = new SolidBrush(Color.FromArgb(r, 0, b)))
                            {
                                g.FillRectangle(brush, legendX, legendY + i, legendWidth, 1);
                            }
                        }
                        g.DrawString("Value", new Font("Consolas", 8F), Brushes.Black, legendX, legendY - 20);
                        g.DrawString("Max", new Font("Consolas", 8F), Brushes.Black, legendX, legendY - 10);
                        g.DrawString("Min", new Font("Consolas", 8F), Brushes.Black, legendX, legendY + legendHeight);
                    }
                    else // Velocity Profile
                    {
                        const int depthLevels = 10;
                        double maxValue = spatialMode.Length > 0 ? Math.Max(Math.Abs(spatialMode.Min()), Math.Abs(spatialMode.Max())) : 1.0;
                        if (maxValue < 1e-10) maxValue = 1.0;

                        g.DrawString($"EOF Mode {modeIndex + 1} for Velocity Profile", new Font("Consolas", 10F, FontStyle.Bold), Brushes.Black, margin, 10);
                        g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin + plotWidth, yOrigin);
                        g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin, yOrigin - plotHeight);
                        g.DrawString("Velocity", new Font("Consolas", 9F), Brushes.Black, xOrigin + plotWidth / 2 - 20, yOrigin + 20);
                        g.DrawString("Depth", new Font("Consolas", 9F), Brushes.Black, xOrigin - 30, yOrigin - plotHeight / 2 - 10, new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical });

                        DrawGridAndTicks(g, xOrigin, yOrigin, plotWidth, plotHeight, maxValue, depthLevels, false, margin, false, 0);

                        for (int i = 1; i < depthLevels; i++)
                        {
                            float x1 = xOrigin + (float)(spatialMode[i - 1] / maxValue * plotWidth / 2 + plotWidth / 2);
                            float y1 = yOrigin - (float)((i - 1) / (double)(depthLevels - 1) * plotHeight);
                            float x2 = xOrigin + (float)(spatialMode[i] / maxValue * plotWidth / 2 + plotWidth / 2);
                            float y2 = yOrigin - (float)(i / (double)(depthLevels - 1) * plotHeight);
                            x1 = Math.Max(margin, Math.Min(xOrigin + plotWidth, x1));
                            x2 = Math.Max(margin, Math.Min(xOrigin + plotWidth, x2));
                            y1 = Math.Max(margin, Math.Min(yOrigin, y1));
                            y2 = Math.Max(margin, Math.Min(yOrigin, y2));
                            g.DrawLine(Pens.Blue, x1, y1, x2, y2);
                        }

                        g.DrawLine(Pens.Blue, xOrigin + plotWidth - 60, 30, xOrigin + plotWidth - 40, 30);
                        g.DrawString("Velocity Profile", new Font("Consolas", 8F), Brushes.Black, xOrigin + plotWidth - 100, 25);
                    }
                }
            }
        }

        private void HeatmapPanel_Paint(object sender, PaintEventArgs e)
        {
            if (IsDisposed || heatmapPanel.IsDisposed) return;

            Graphics g = e.Graphics;
            g.Clear(Color.White);
            g.SmoothingMode = System.Drawing.Drawing2D.SmoothingMode.AntiAlias;

            const int margin = 50;
            float plotWidth = heatmapPanel.Width - 2 * margin;
            float plotHeight = heatmapPanel.Height - 2 * margin - 30;
            float xOrigin = margin;
            float yOrigin = heatmapPanel.Height - margin;

            bool isLogLog = scaleComboBox.SelectedItem?.ToString() == "Log-Log";
            string analysisType = analysisComboBox.SelectedItem?.ToString() ?? "Spectral Analysis";

            if (analysisType == "Spectral Analysis")
            {
                if (heatmapPanel.Tag is HeatmapData heatmapData)
                {
                    double[,] powerSpectrum = heatmapData.PowerSpectrum;
                    double[] frequencies = heatmapData.Frequencies;
                    double[] times = heatmapData.Times;
                    int segmentLength = heatmapData.SegmentLength;
                    double maxPower = 0;
                    for (int i = 1; i < frequencies.Length; i++)
                        for (int j = 0; j < times.Length; j++)
                            if (powerSpectrum[i, j] > maxPower) maxPower = powerSpectrum[i, j];
                    double maxFreq = frequencies.Length > 0 ? frequencies[frequencies.Length - 1] : 1.0;
                    double minFreq = frequencies.Length > 0 ? frequencies[1] : 1e-6;
                    double maxTime = times.Length > 0 ? times[times.Length - 1] : duration;

                    if (maxPower < 1e-10) maxPower = 1.0;
                    if (minFreq < 1e-6) minFreq = 1e-6;
                    if (maxFreq <= minFreq) maxFreq = minFreq + 1e-6;
                    if (maxTime <= 0) maxTime = duration;

                    g.DrawString($"Time-Frequency Heatmap of {variableComboBox.SelectedItem} ({(isLogLog ? "Log-Log" : "Linear")}, {selectedWindowType} Window)", new Font("Consolas", 10F, FontStyle.Bold), Brushes.Black, margin, 10);
                    g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin + plotWidth, yOrigin);
                    g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin, yOrigin - plotHeight);
                    g.DrawString("Time (s)", new Font("Consolas", 9F), Brushes.Black, xOrigin + plotWidth / 2 - 20, yOrigin + 20);
                    g.DrawString("Frequency (Hz)", new Font("Consolas", 9F), Brushes.Black, xOrigin - 40, yOrigin - plotHeight / 2 - 10, new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical });

                    DrawGridAndTicks(g, xOrigin, yOrigin, plotWidth, plotHeight, maxTime, maxFreq, true, margin, isLogLog, minFreq);

                    float cellWidth = plotWidth / times.Length;
                    float cellHeight = plotHeight / (frequencies.Length - 1);
                    for (int j = 0; j < times.Length; j++)
                    {
                        for (int i = 1; i < frequencies.Length; i++)
                        {
                            if (frequencies[i] <= 0 || double.IsNaN(powerSpectrum[i, j])) continue;

                            double power = isLogLog ? Math.Log10(Math.Max(powerSpectrum[i, j], 1e-10)) : powerSpectrum[i, j];
                            double maxPowerScaled = isLogLog ? Math.Log10(maxPower) : maxPower;
                            double minPowerScaled = isLogLog ? Math.Log10(1e-10) : 0;
                            double freqScale = isLogLog ? (Math.Log10(frequencies[i]) - Math.Log10(minFreq)) / (Math.Log10(maxFreq) - Math.Log10(minFreq)) : frequencies[i] / maxFreq;
                            double powerScale = maxPowerScaled > minPowerScaled ? (power - minPowerScaled) / (maxPowerScaled - minPowerScaled) : 0;

                            if (double.IsNaN(powerScale) || double.IsInfinity(powerScale)) powerScale = 0;
                            powerScale = Math.Max(0, Math.Min(1, powerScale));

                            float x = xOrigin + j * cellWidth;
                            float y = yOrigin - (float)(freqScale * plotHeight);
                            float width = cellWidth;
                            float height = cellHeight;
                            if (isLogLog)
                            {
                                height = plotHeight / (frequencies.Length - 1);
                            }

                            x = Math.Max(margin, Math.Min(xOrigin + plotWidth, x));
                            y = Math.Max(margin, Math.Min(yOrigin, y));
                            width = Math.Max(1, Math.Min(plotWidth - (x - xOrigin), width));
                            height = Math.Max(1, Math.Min(yOrigin - y, height));

                            int r = (int)(powerScale * 255);
                            int b = (int)((1 - powerScale) * 255);
                            using (SolidBrush brush = new SolidBrush(Color.FromArgb(r, 0, b)))
                            {
                                g.FillRectangle(brush, x, y, width, height);
                            }
                        }
                    }

                    double m2Freq = 1.0 / 44712;
                    double k1Freq = 1.0 / 86148;
                    double freqTolerance = (1.0 / (samplingRate * segmentLength)) * 2;
                    foreach (var freq in new[] { (m2Freq, "M2"), (k1Freq, "K1") })
                    {
                        if (freq.Item1 <= maxFreq && freq.Item1 >= minFreq)
                        {
                            double yScale = isLogLog ? (Math.Log10(freq.Item1) - Math.Log10(minFreq)) / (Math.Log10(maxFreq) - Math.Log10(minFreq)) : freq.Item1 / maxFreq;
                            if (double.IsNaN(yScale) || double.IsInfinity(yScale))
                                continue;

                            float yPos = yOrigin - (float)(yScale * plotHeight);
                            yPos = Math.Max(margin, Math.Min(yOrigin, yPos));
                            g.DrawLine(Pens.Red, xOrigin, yPos, xOrigin + plotWidth, yPos);
                            g.DrawString(freq.Item2, new Font("Consolas", 8F), Brushes.Red, xOrigin + plotWidth - 20, yPos - 10);
                        }
                    }

                    float legendWidth = 20;
                    float legendHeight = plotHeight / 2;
                    float legendX = xOrigin + plotWidth - 30;
                    float legendY = yOrigin - legendHeight;
                    for (int i = 0; i < (int)legendHeight; i++)
                    {
                        float powerScale = i / legendHeight;
                        int r = (int)(powerScale * 255);
                        int b = (int)((1 - powerScale) * 255);
                        using (SolidBrush brush = new SolidBrush(Color.FromArgb(r, 0, b)))
                        {
                            g.FillRectangle(brush, legendX, legendY + i, legendWidth, 1);
                        }
                    }
                    g.DrawString("Power", new Font("Consolas", 8F), Brushes.Black, legendX, legendY - 20);
                    g.DrawString(isLogLog ? "Log" : "Max", new Font("Consolas", 8F), Brushes.Black, legendX, legendY - 10);
                    g.DrawString("Min", new Font("Consolas", 8F), Brushes.Black, legendX, legendY + legendHeight);
                }
            }
            else // EOF/POD Analysis
            {
                if (heatmapPanel.Tag is EOFData eofData)
                {
                    int modeIndex = modeComboBox.SelectedIndex;
                    if (modeIndex < 0 || modeIndex >= eofData.NumModes) modeIndex = 0;
                    double[] temporalCoefficients = eofData.TemporalCoefficients[modeIndex];
                    double maxValue = temporalCoefficients.Length > 0 ? Math.Max(Math.Abs(temporalCoefficients.Min()), Math.Abs(temporalCoefficients.Max())) : 1.0;
                    double maxTime = duration;

                    if (maxValue < 1e-10) maxValue = 1.0;

                    g.DrawString($"Temporal Coefficient for Mode {modeIndex + 1} ({variableComboBox.SelectedItem})", new Font("Consolas", 10F, FontStyle.Bold), Brushes.Black, margin, 10);
                    g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin + plotWidth, yOrigin);
                    g.DrawLine(Pens.Black, xOrigin, yOrigin, xOrigin, yOrigin - plotHeight);
                    g.DrawString("Time (s)", new Font("Consolas", 9F), Brushes.Black, xOrigin + plotWidth / 2 - 20, yOrigin + 20);
                    g.DrawString("Coefficient", new Font("Consolas", 9F), Brushes.Black, xOrigin - 40, yOrigin - plotHeight / 2 - 10, new StringFormat { FormatFlags = StringFormatFlags.DirectionVertical });

                    DrawGridAndTicks(g, xOrigin, yOrigin, plotWidth, plotHeight, maxTime, maxValue, false, margin, false, 0);

                    for (int i = 1; i < temporalCoefficients.Length; i++)
                    {
                        float x1 = xOrigin + (float)((i - 1) * samplingRate / maxTime * plotWidth);
                        float y1 = yOrigin - (float)(temporalCoefficients[i - 1] / maxValue * plotHeight / 2 + plotHeight / 2);
                        float x2 = xOrigin + (float)(i * samplingRate / maxTime * plotWidth);
                        float y2 = yOrigin - (float)(temporalCoefficients[i] / maxValue * plotHeight / 2 + plotHeight / 2);
                        y1 = Math.Max(margin, Math.Min(yOrigin, y1));
                        y2 = Math.Max(margin, Math.Min(yOrigin, y2));
                        g.DrawLine(Pens.Blue, x1, y1, x2, y2);
                    }

                    g.DrawLine(Pens.Blue, xOrigin + plotWidth - 60, 30, xOrigin + plotWidth - 40, 30);
                    g.DrawString("Temporal Coefficient", new Font("Consolas", 8F), Brushes.Black, xOrigin + plotWidth - 100, 25);
                }
            }
        }

        private void DrawGridAndTicks(Graphics g, float xOrigin, float yOrigin, float plotWidth, float plotHeight, double maxX, double maxY, bool isPowerSpectrum, int margin, bool isLogLog, double minFreq)
        {
            if (maxX < 1E-10) maxX = 1.0;
            if (maxY < 1E-10) maxY = 1.0;
            if (isPowerSpectrum && isLogLog && minFreq < 1e-6) minFreq = 1e-6;

            double xTickInterval = isPowerSpectrum && isLogLog ? CalculateLogTickInterval(minFreq, maxX, 5) : CalculateTickInterval(maxX, 5);
            double yTickInterval = isPowerSpectrum && isLogLog ? CalculateLogTickInterval(1e-10, maxY, 5) : CalculateTickInterval(maxY, 5);

            // Draw X-axis grid and ticks
            if (isPowerSpectrum && isLogLog)
            {
                double logMinFreq = Math.Log10(minFreq);
                double logMaxFreq = Math.Log10(maxX);
                for (double logX = Math.Ceiling(logMinFreq); logX <= logMaxFreq; logX += xTickInterval)
                {
                    double x = Math.Pow(10, logX);
                    float xPos = xOrigin + (float)((logX - logMinFreq) / (logMaxFreq - logMinFreq) * plotWidth);
                    xPos = Math.Max(margin, Math.Min(xOrigin + plotWidth, xPos));
                    g.DrawLine(Pens.LightGray, xPos, yOrigin, xPos, yOrigin - plotHeight);
                    g.DrawLine(Pens.Black, xPos, yOrigin, xPos, yOrigin + 5);
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
                    g.DrawLine(Pens.LightGray, xPos, yOrigin, xPos, yOrigin - plotHeight);
                    g.DrawLine(Pens.Black, xPos, yOrigin, xPos, yOrigin + 5);
                    string label = isPowerSpectrum ? x.ToString("F4") : ((int)x).ToString();
                    g.DrawString(label, new Font("Consolas", 8F), Brushes.Black, xPos - 15, yOrigin + 10);
                }
            }

            // Draw Y-axis grid and ticks
            if (isPowerSpectrum && isLogLog)
            {
                double logMinY = Math.Log10(1e-10); // Minimum power for log scale
                double logMaxY = Math.Log10(maxY);
                for (double logY = Math.Ceiling(logMinY); logY <= logMaxY; logY += yTickInterval)
                {
                    double y = Math.Pow(10, logY);
                    float yPos = yOrigin - (float)((logY - logMinY) / (logMaxY - logMinY) * plotHeight);
                    yPos = Math.Max(margin, Math.Min(yOrigin, yPos));
                    g.DrawLine(Pens.LightGray, xOrigin, yPos, xOrigin + plotWidth, yPos);
                    g.DrawLine(Pens.Black, xOrigin - 5, yPos, xOrigin, yPos);
                    string label = $"10^{logY:F0}";
                    g.DrawString(label, new Font("Consolas", 8F), Brushes.Black, xOrigin - 40, yPos - 5);
                }
            }
            else
            {
                double yStart = isPowerSpectrum ? 0 : -maxY;
                double yEnd = maxY;
                for (double y = yStart; y <= yEnd; y += yTickInterval)
                {
                    float yPos = isPowerSpectrum
                        ? yOrigin - (float)(y / maxY * plotHeight)
                        : yOrigin - (float)((y + maxY) / (2 * maxY) * plotHeight);
                    yPos = Math.Max(margin, Math.Min(yOrigin, yPos));
                    g.DrawLine(Pens.LightGray, xOrigin, yPos, xOrigin + plotWidth, yPos);
                    g.DrawLine(Pens.Black, xOrigin - 5, yPos, xOrigin, yPos);
                    string label = isPowerSpectrum ? y.ToString("F2") : y.ToString("F1");
                    g.DrawString(label, new Font("Consolas", 8F), Brushes.Black, xOrigin - 40, yPos - 5);
                }
            }
        }
        private void FourierTransform(Complex[] data)
        {
            int n = data.Length;
            if (n <= 1) return;

            // Bit-reversal permutation
            int logN = (int)Math.Log(n, 2);
            if (Math.Pow(2, logN) != n) throw new ArgumentException("Data length must be a power of 2.");

            for (int i = 0; i < n; i++)
            {
                int rev = 0;
                for (int j = 0; j < logN; j++)
                {
                    if ((i & (1 << j)) != 0) rev |= (1 << (logN - 1 - j));
                }
                if (rev > i)
                {
                    Complex temp = data[i];
                    data[i] = data[rev];
                    data[rev] = temp;
                }
            }

            // Cooley-Tukey FFT
            for (int size = 2; size <= n; size *= 2)
            {
                int halfSize = size / 2;
                double angleStep = -2 * Math.PI / size;
                for (int i = 0; i < n; i += size)
                {
                    for (int j = i; j < i + halfSize; j++)
                    {
                        int k = j + halfSize;
                        Complex w = new Complex(Math.Cos(angleStep * (j - i)), Math.Sin(angleStep * (j - i)));
                        Complex t = w * data[k];
                        data[k] = data[j] - t;
                        data[j] = data[j] + t;
                    }
                }
            }
        }

        private double CalculateTickInterval(double range, int numTicks)
        {
            if (range <= 1e-10 || numTicks <= 0) return 1.0;
            double roughInterval = range / numTicks;
            double magnitude = Math.Pow(10, Math.Floor(Math.Log10(roughInterval)));
            double normalized = roughInterval / magnitude;
            double step = normalized <= 1.5 ? 1 : normalized <= 3 ? 2 : normalized <= 7 ? 5 : 10;
            return step * magnitude;
        }

        private double CalculateLogTickInterval(double min, double max, int numTicks)
        {
            if (min <= 0 || max <= min || numTicks <= 0) return 1.0;
            double logMin = Math.Log10(min);
            double logMax = Math.Log10(max);
            return (logMax - logMin) / numTicks;
        }

        protected override void Dispose(bool disposing)
        {
            if (disposing)
            {
                simulationTimer?.Dispose();
                visualizationPanel?.Dispose();
                heatmapPanel?.Dispose();
                outputTextBox?.Dispose();
                variableComboBox?.Dispose();
                scaleComboBox?.Dispose();
                windowComboBox?.Dispose();
                analysisComboBox?.Dispose();
                modeComboBox?.Dispose();
                analyzeButton?.Dispose();
                pauseButton?.Dispose();
                resetButton?.Dispose();
                clearButton?.Dispose();
            }
            base.Dispose(disposing);
        }
    }
}
