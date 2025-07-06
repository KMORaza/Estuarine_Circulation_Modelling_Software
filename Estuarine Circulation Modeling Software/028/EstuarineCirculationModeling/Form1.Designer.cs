using System.Windows.Forms;

namespace EstuarineCirculationModeling
{
    partial class Form1
    {
        private System.ComponentModel.IContainer components = null;

        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        private void InitializeComponent()
        {
            this.controlPanel = new System.Windows.Forms.Panel();
            this.tidalStrainButton = new System.Windows.Forms.Button();
            this.satButton = new System.Windows.Forms.Button();
            this.compForcingButton = new System.Windows.Forms.Button();
            this.tamButton = new System.Windows.Forms.Button();
            this.atmButton = new System.Windows.Forms.Button(); // Added for ATM
            this.lbLesButton = new System.Windows.Forms.Button();
            this.largeEddyButton = new System.Windows.Forms.Button();
            this.stratificationButton = new System.Windows.Forms.Button();
            this.useRANSSolverCheckBox = new System.Windows.Forms.CheckBox();
            this.useNonHydrostaticCheckBox = new System.Windows.Forms.CheckBox();
            this.startButton = new System.Windows.Forms.Button();
            this.pauseButton = new System.Windows.Forms.Button();
            this.resetButton = new System.Windows.Forms.Button();
            this.riverInflowLabel = new System.Windows.Forms.Label();
            this.riverInflowTextBox = new System.Windows.Forms.TextBox();
            this.tidalAmplitudeLabel = new System.Windows.Forms.Label();
            this.tidalAmplitudeTextBox = new System.Windows.Forms.TextBox();
            this.tidalPeriodLabel = new System.Windows.Forms.Label();
            this.tidalPeriodTextBox = new System.Windows.Forms.TextBox();
            this.oceanSalinityLabel = new System.Windows.Forms.Label();
            this.oceanSalinityTextBox = new System.Windows.Forms.TextBox();
            this.oceanTemperatureLabel = new System.Windows.Forms.Label();
            this.oceanTemperatureTextBox = new System.Windows.Forms.TextBox();
            this.estuaryLengthLabel = new System.Windows.Forms.Label();
            this.estuaryLengthTextBox = new System.Windows.Forms.TextBox();
            this.estuaryDepthLabel = new System.Windows.Forms.Label();
            this.estuaryDepthTextBox = new System.Windows.Forms.TextBox();
            this.outputConsoleTextBox = new System.Windows.Forms.TextBox();
            this.visualizationPanel = new System.Windows.Forms.Panel();
            this.controlPanel.SuspendLayout();
            this.SuspendLayout();

            // controlPanel
            this.controlPanel.Controls.Add(this.tidalStrainButton);
            this.controlPanel.Controls.Add(this.satButton);
            this.controlPanel.Controls.Add(this.compForcingButton);
            this.controlPanel.Controls.Add(this.tamButton);
            this.controlPanel.Controls.Add(this.atmButton); // Added for ATM
            this.controlPanel.Controls.Add(this.lbLesButton);
            this.controlPanel.Controls.Add(this.largeEddyButton);
            this.controlPanel.Controls.Add(this.stratificationButton);
            this.controlPanel.Controls.Add(this.useRANSSolverCheckBox);
            this.controlPanel.Controls.Add(this.useNonHydrostaticCheckBox);
            this.controlPanel.Controls.Add(this.startButton);
            this.controlPanel.Controls.Add(this.pauseButton);
            this.controlPanel.Controls.Add(this.resetButton);
            this.controlPanel.Controls.Add(this.riverInflowLabel);
            this.controlPanel.Controls.Add(this.riverInflowTextBox);
            this.controlPanel.Controls.Add(this.tidalAmplitudeLabel);
            this.controlPanel.Controls.Add(this.tidalAmplitudeTextBox);
            this.controlPanel.Controls.Add(this.tidalPeriodLabel);
            this.controlPanel.Controls.Add(this.tidalPeriodTextBox);
            this.controlPanel.Controls.Add(this.oceanSalinityLabel);
            this.controlPanel.Controls.Add(this.oceanSalinityTextBox);
            this.controlPanel.Controls.Add(this.oceanTemperatureLabel);
            this.controlPanel.Controls.Add(this.oceanTemperatureTextBox);
            this.controlPanel.Controls.Add(this.estuaryLengthLabel);
            this.controlPanel.Controls.Add(this.estuaryLengthTextBox);
            this.controlPanel.Controls.Add(this.estuaryDepthLabel);
            this.controlPanel.Controls.Add(this.estuaryDepthTextBox);
            this.controlPanel.Location = new System.Drawing.Point(12, 12);
            this.controlPanel.Name = "controlPanel";
            this.controlPanel.Size = new System.Drawing.Size(300, 520);
            this.controlPanel.TabIndex = 0;
            this.controlPanel.AutoScroll = true;
            this.controlPanel.AutoSize = false;
            this.controlPanel.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.controlPanel.Font = new System.Drawing.Font("Consolas", 9F);
            this.controlPanel.BackColor = System.Drawing.SystemColors.Control;
            this.controlPanel.ForeColor = System.Drawing.Color.Black;

            // riverInflowLabel
            this.riverInflowLabel.AutoSize = true;
            this.riverInflowLabel.Location = new System.Drawing.Point(10, 20);
            this.riverInflowLabel.Name = "riverInflowLabel";
            this.riverInflowLabel.Size = new System.Drawing.Size(80, 15);
            this.riverInflowLabel.TabIndex = 0;
            this.riverInflowLabel.Text = "River Inflow (m³/s):";
            this.riverInflowLabel.Font = new System.Drawing.Font("Consolas", 9F);
            this.riverInflowLabel.BackColor = System.Drawing.SystemColors.Control;
            this.riverInflowLabel.ForeColor = System.Drawing.Color.Black;

            // riverInflowTextBox
            this.riverInflowTextBox.Location = new System.Drawing.Point(10, 40);
            this.riverInflowTextBox.Name = "riverInflowTextBox";
            this.riverInflowTextBox.Size = new System.Drawing.Size(250, 22);
            this.riverInflowTextBox.TabIndex = 1;
            this.riverInflowTextBox.Text = "0.1";
            this.riverInflowTextBox.Font = new System.Drawing.Font("Consolas", 9F);
            this.riverInflowTextBox.BackColor = System.Drawing.Color.White;
            this.riverInflowTextBox.ForeColor = System.Drawing.Color.Black;
            this.riverInflowTextBox.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;

            // tidalAmplitudeLabel
            this.tidalAmplitudeLabel.AutoSize = true;
            this.tidalAmplitudeLabel.Location = new System.Drawing.Point(10, 70);
            this.tidalAmplitudeLabel.Name = "tidalAmplitudeLabel";
            this.tidalAmplitudeLabel.Size = new System.Drawing.Size(80, 15);
            this.tidalAmplitudeLabel.TabIndex = 2;
            this.tidalAmplitudeLabel.Text = "Tidal Amplitude (m):";
            this.tidalAmplitudeLabel.Font = new System.Drawing.Font("Consolas", 9F);
            this.tidalAmplitudeLabel.BackColor = System.Drawing.SystemColors.Control;
            this.tidalAmplitudeLabel.ForeColor = System.Drawing.Color.Black;

            // tidalAmplitudeTextBox
            this.tidalAmplitudeTextBox.Location = new System.Drawing.Point(10, 90);
            this.tidalAmplitudeTextBox.Name = "tidalAmplitudeTextBox";
            this.tidalAmplitudeTextBox.Size = new System.Drawing.Size(250, 22);
            this.tidalAmplitudeTextBox.TabIndex = 3;
            this.tidalAmplitudeTextBox.Text = "1.0";
            this.tidalAmplitudeTextBox.Font = new System.Drawing.Font("Consolas", 9F);
            this.tidalAmplitudeTextBox.BackColor = System.Drawing.Color.White;
            this.tidalAmplitudeTextBox.ForeColor = System.Drawing.Color.Black;
            this.tidalAmplitudeTextBox.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;

            // tidalPeriodLabel
            this.tidalPeriodLabel.AutoSize = true;
            this.tidalPeriodLabel.Location = new System.Drawing.Point(10, 120);
            this.tidalPeriodLabel.Name = "tidalPeriodLabel";
            this.tidalPeriodLabel.Size = new System.Drawing.Size(80, 15);
            this.tidalPeriodLabel.TabIndex = 4;
            this.tidalPeriodLabel.Text = "Tidal Period (s):";
            this.tidalPeriodLabel.Font = new System.Drawing.Font("Consolas", 9F);
            this.tidalPeriodLabel.BackColor = System.Drawing.SystemColors.Control;
            this.tidalPeriodLabel.ForeColor = System.Drawing.Color.Black;

            // tidalPeriodTextBox
            this.tidalPeriodTextBox.Location = new System.Drawing.Point(10, 140);
            this.tidalPeriodTextBox.Name = "tidalPeriodTextBox";
            this.tidalPeriodTextBox.Size = new System.Drawing.Size(250, 22);
            this.tidalPeriodTextBox.TabIndex = 5;
            this.tidalPeriodTextBox.Text = "43200";
            this.tidalPeriodTextBox.Font = new System.Drawing.Font("Consolas", 9F);
            this.tidalPeriodTextBox.BackColor = System.Drawing.Color.White;
            this.tidalPeriodTextBox.ForeColor = System.Drawing.Color.Black;
            this.tidalPeriodTextBox.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;

            // oceanSalinityLabel
            this.oceanSalinityLabel.AutoSize = true;
            this.oceanSalinityLabel.Location = new System.Drawing.Point(10, 170);
            this.oceanSalinityLabel.Name = "oceanSalinityLabel";
            this.oceanSalinityLabel.Size = new System.Drawing.Size(80, 15);
            this.oceanSalinityLabel.TabIndex = 6;
            this.oceanSalinityLabel.Text = "Ocean Salinity (PSU):";
            this.oceanSalinityLabel.Font = new System.Drawing.Font("Consolas", 9F);
            this.oceanSalinityLabel.BackColor = System.Drawing.SystemColors.Control;
            this.oceanSalinityLabel.ForeColor = System.Drawing.Color.Black;

            // oceanSalinityTextBox
            this.oceanSalinityTextBox.Location = new System.Drawing.Point(10, 190);
            this.oceanSalinityTextBox.Name = "oceanSalinityTextBox";
            this.oceanSalinityTextBox.Size = new System.Drawing.Size(250, 22);
            this.oceanSalinityTextBox.TabIndex = 7;
            this.oceanSalinityTextBox.Text = "35.0";
            this.oceanSalinityTextBox.Font = new System.Drawing.Font("Consolas", 9F);
            this.oceanSalinityTextBox.BackColor = System.Drawing.Color.White;
            this.oceanSalinityTextBox.ForeColor = System.Drawing.Color.Black;
            this.oceanSalinityTextBox.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;

            // oceanTemperatureLabel
            this.oceanTemperatureLabel.AutoSize = true;
            this.oceanTemperatureLabel.Location = new System.Drawing.Point(10, 220);
            this.oceanTemperatureLabel.Name = "oceanTemperatureLabel";
            this.oceanTemperatureLabel.Size = new System.Drawing.Size(100, 15);
            this.oceanTemperatureLabel.TabIndex = 8;
            this.oceanTemperatureLabel.Text = "Ocean Temperature (°C):";
            this.oceanTemperatureLabel.Font = new System.Drawing.Font("Consolas", 9F);
            this.oceanTemperatureLabel.BackColor = System.Drawing.SystemColors.Control;
            this.oceanTemperatureLabel.ForeColor = System.Drawing.Color.Black;

            // oceanTemperatureTextBox
            this.oceanTemperatureTextBox.Location = new System.Drawing.Point(10, 240);
            this.oceanTemperatureTextBox.Name = "oceanTemperatureTextBox";
            this.oceanTemperatureTextBox.Size = new System.Drawing.Size(250, 22);
            this.oceanTemperatureTextBox.TabIndex = 9;
            this.oceanTemperatureTextBox.Text = "20.0";
            this.oceanTemperatureTextBox.Font = new System.Drawing.Font("Consolas", 9F);
            this.oceanTemperatureTextBox.BackColor = System.Drawing.Color.White;
            this.oceanTemperatureTextBox.ForeColor = System.Drawing.Color.Black;
            this.oceanTemperatureTextBox.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;

            // estuaryLengthLabel
            this.estuaryLengthLabel.AutoSize = true;
            this.estuaryLengthLabel.Location = new System.Drawing.Point(10, 270);
            this.estuaryLengthLabel.Name = "estuaryLengthLabel";
            this.estuaryLengthLabel.Size = new System.Drawing.Size(80, 15);
            this.estuaryLengthLabel.TabIndex = 10;
            this.estuaryLengthLabel.Text = "Estuary Length (m):";
            this.estuaryLengthLabel.Font = new System.Drawing.Font("Consolas", 9F);
            this.estuaryLengthLabel.BackColor = System.Drawing.SystemColors.Control;
            this.estuaryLengthLabel.ForeColor = System.Drawing.Color.Black;

            // estuaryLengthTextBox
            this.estuaryLengthTextBox.Location = new System.Drawing.Point(10, 290);
            this.estuaryLengthTextBox.Name = "estuaryLengthTextBox";
            this.estuaryLengthTextBox.Size = new System.Drawing.Size(250, 22);
            this.estuaryLengthTextBox.TabIndex = 11;
            this.estuaryLengthTextBox.Text = "10000";
            this.estuaryLengthTextBox.Font = new System.Drawing.Font("Consolas", 9F);
            this.estuaryLengthTextBox.BackColor = System.Drawing.Color.White;
            this.estuaryLengthTextBox.ForeColor = System.Drawing.Color.Black;
            this.estuaryLengthTextBox.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;

            // estuaryDepthLabel
            this.estuaryDepthLabel.AutoSize = true;
            this.estuaryDepthLabel.Location = new System.Drawing.Point(10, 320);
            this.estuaryDepthLabel.Name = "estuaryDepthLabel";
            this.estuaryDepthLabel.Size = new System.Drawing.Size(80, 15);
            this.estuaryDepthLabel.TabIndex = 12;
            this.estuaryDepthLabel.Text = "Estuary Depth (m):";
            this.estuaryDepthLabel.Font = new System.Drawing.Font("Consolas", 9F);
            this.estuaryDepthLabel.BackColor = System.Drawing.SystemColors.Control;
            this.estuaryDepthLabel.ForeColor = System.Drawing.Color.Black;

            // estuaryDepthTextBox
            this.estuaryDepthTextBox.Location = new System.Drawing.Point(10, 340);
            this.estuaryDepthTextBox.Name = "estuaryDepthTextBox";
            this.estuaryDepthTextBox.Size = new System.Drawing.Size(250, 22);
            this.estuaryDepthTextBox.TabIndex = 13;
            this.estuaryDepthTextBox.Text = "10.0";
            this.estuaryDepthTextBox.Font = new System.Drawing.Font("Consolas", 9F);
            this.estuaryDepthTextBox.BackColor = System.Drawing.Color.White;
            this.estuaryDepthTextBox.ForeColor = System.Drawing.Color.Black;
            this.estuaryDepthTextBox.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;

            // useRANSSolverCheckBox
            this.useRANSSolverCheckBox.AutoSize = true;
            this.useRANSSolverCheckBox.Location = new System.Drawing.Point(10, 370);
            this.useRANSSolverCheckBox.Name = "useRANSSolverCheckBox";
            this.useRANSSolverCheckBox.Size = new System.Drawing.Size(120, 19);
            this.useRANSSolverCheckBox.TabIndex = 14;
            this.useRANSSolverCheckBox.Text = "Use RANS Solver";
            this.useRANSSolverCheckBox.Font = new System.Drawing.Font("Consolas", 9F);
            this.useRANSSolverCheckBox.BackColor = System.Drawing.SystemColors.Control;
            this.useRANSSolverCheckBox.ForeColor = System.Drawing.Color.Black;
            this.useRANSSolverCheckBox.CheckedChanged += new System.EventHandler(this.useRANSSolverCheckBox_CheckedChanged);

            // useNonHydrostaticCheckBox
            this.useNonHydrostaticCheckBox.AutoSize = true;
            this.useNonHydrostaticCheckBox.Location = new System.Drawing.Point(10, 390);
            this.useNonHydrostaticCheckBox.Name = "useNonHydrostaticCheckBox";
            this.useNonHydrostaticCheckBox.Size = new System.Drawing.Size(140, 19);
            this.useNonHydrostaticCheckBox.TabIndex = 15;
            this.useNonHydrostaticCheckBox.Text = "Use Non-Hydrostatic";
            this.useNonHydrostaticCheckBox.Enabled = false;
            this.useNonHydrostaticCheckBox.Font = new System.Drawing.Font("Consolas", 9F);
            this.useNonHydrostaticCheckBox.BackColor = System.Drawing.SystemColors.Control;
            this.useNonHydrostaticCheckBox.ForeColor = System.Drawing.Color.Black;
            this.useNonHydrostaticCheckBox.CheckedChanged += new System.EventHandler(this.useNonHydrostaticCheckBox_CheckedChanged);

            // stratificationButton
            this.stratificationButton.Location = new System.Drawing.Point(10, 420);
            this.stratificationButton.Name = "stratificationButton";
            this.stratificationButton.Size = new System.Drawing.Size(170, 25);
            this.stratificationButton.TabIndex = 16;
            this.stratificationButton.Text = "Stratification";
            this.stratificationButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.stratificationButton.FlatAppearance.BorderSize = 1;
            this.stratificationButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.stratificationButton.BackColor = System.Drawing.SystemColors.Control;
            this.stratificationButton.ForeColor = System.Drawing.Color.Black;
            this.stratificationButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.stratificationButton.Click += new System.EventHandler(this.stratificationButton_Click);

            // largeEddyButton
            this.largeEddyButton.Location = new System.Drawing.Point(10, 455);
            this.largeEddyButton.Name = "largeEddyButton";
            this.largeEddyButton.Size = new System.Drawing.Size(170, 25);
            this.largeEddyButton.TabIndex = 17;
            this.largeEddyButton.Text = "Large Eddy";
            this.largeEddyButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.largeEddyButton.FlatAppearance.BorderSize = 1;
            this.largeEddyButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.largeEddyButton.BackColor = System.Drawing.SystemColors.Control;
            this.largeEddyButton.ForeColor = System.Drawing.Color.Black;
            this.largeEddyButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.largeEddyButton.Click += new System.EventHandler(this.largeEddyButton_Click);

            // lbLesButton
            this.lbLesButton.Location = new System.Drawing.Point(10, 490);
            this.lbLesButton.Name = "lbLesButton";
            this.lbLesButton.Size = new System.Drawing.Size(170, 25);
            this.lbLesButton.TabIndex = 18;
            this.lbLesButton.Text = "LB-LES";
            this.lbLesButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.lbLesButton.FlatAppearance.BorderSize = 1;
            this.lbLesButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.lbLesButton.BackColor = System.Drawing.SystemColors.Control;
            this.lbLesButton.ForeColor = System.Drawing.Color.Black;
            this.lbLesButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.lbLesButton.Click += new System.EventHandler(this.lbLesButton_Click);

            // tidalStrainButton
            this.tidalStrainButton.Location = new System.Drawing.Point(10, 525);
            this.tidalStrainButton.Name = "tidalStrainButton";
            this.tidalStrainButton.Size = new System.Drawing.Size(170, 25);
            this.tidalStrainButton.TabIndex = 24;
            this.tidalStrainButton.Text = "Tidal Strain";
            this.tidalStrainButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.tidalStrainButton.FlatAppearance.BorderSize = 1;
            this.tidalStrainButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.tidalStrainButton.BackColor = System.Drawing.SystemColors.Control;
            this.tidalStrainButton.ForeColor = System.Drawing.Color.Black;
            this.tidalStrainButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.tidalStrainButton.Click += new System.EventHandler(this.tidalStrainButton_Click);

            // satButton
            this.satButton.Location = new System.Drawing.Point(10, 560);
            this.satButton.Name = "satButton";
            this.satButton.Size = new System.Drawing.Size(170, 25);
            this.satButton.TabIndex = 25;
            this.satButton.Text = "SAT";
            this.satButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.satButton.FlatAppearance.BorderSize = 1;
            this.satButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.satButton.BackColor = System.Drawing.SystemColors.Control;
            this.satButton.ForeColor = System.Drawing.Color.Black;
            this.satButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.satButton.Click += new System.EventHandler(this.satButton_Click);

            // compForcingButton
            this.compForcingButton.Location = new System.Drawing.Point(10, 595);
            this.compForcingButton.Name = "compForcingButton";
            this.compForcingButton.Size = new System.Drawing.Size(170, 25);
            this.compForcingButton.TabIndex = 26;
            this.compForcingButton.Text = "Comp. Forcing";
            this.compForcingButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.compForcingButton.FlatAppearance.BorderSize = 1;
            this.compForcingButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.compForcingButton.BackColor = System.Drawing.SystemColors.Control;
            this.compForcingButton.ForeColor = System.Drawing.Color.Black;
            this.compForcingButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.compForcingButton.Click += new System.EventHandler(this.compForcingButton_Click);

            // tamButton
            this.tamButton.Location = new System.Drawing.Point(10, 630);
            this.tamButton.Name = "tamButton";
            this.tamButton.Size = new System.Drawing.Size(170, 25);
            this.tamButton.TabIndex = 27;
            this.tamButton.Text = "Simulate TAM";
            this.tamButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.tamButton.FlatAppearance.BorderSize = 1;
            this.tamButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.tamButton.BackColor = System.Drawing.SystemColors.Control;
            this.tamButton.ForeColor = System.Drawing.Color.Black;
            this.tamButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.tamButton.Click += new System.EventHandler(this.tamButton_Click);

            // atmButton
            this.atmButton.Location = new System.Drawing.Point(10, 665);
            this.atmButton.Name = "atmButton";
            this.atmButton.Size = new System.Drawing.Size(170, 25);
            this.atmButton.TabIndex = 28;
            this.atmButton.Text = "Simulate ATM";
            this.atmButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.atmButton.FlatAppearance.BorderSize = 1;
            this.atmButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.atmButton.BackColor = System.Drawing.SystemColors.Control;
            this.atmButton.ForeColor = System.Drawing.Color.Black;
            this.atmButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.atmButton.Click += new System.EventHandler(this.atmButton_Click);

            // startButton
            this.startButton.Location = new System.Drawing.Point(10, 700);
            this.startButton.Name = "startButton";
            this.startButton.Size = new System.Drawing.Size(170, 25);
            this.startButton.TabIndex = 19;
            this.startButton.Text = "Start";
            this.startButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.startButton.FlatAppearance.BorderSize = 1;
            this.startButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.startButton.BackColor = System.Drawing.SystemColors.Control;
            this.startButton.ForeColor = System.Drawing.Color.Black;
            this.startButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.startButton.Click += new System.EventHandler(this.startButton_Click);

            // pauseButton
            this.pauseButton.Location = new System.Drawing.Point(10, 735);
            this.pauseButton.Name = "pauseButton";
            this.pauseButton.Size = new System.Drawing.Size(170, 25);
            this.pauseButton.TabIndex = 20;
            this.pauseButton.Text = "Pause";
            this.pauseButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.pauseButton.FlatAppearance.BorderSize = 1;
            this.pauseButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.pauseButton.BackColor = System.Drawing.SystemColors.Control;
            this.pauseButton.ForeColor = System.Drawing.Color.Black;
            this.pauseButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.pauseButton.Click += new System.EventHandler(this.pauseButton_Click);

            // resetButton
            this.resetButton.Location = new System.Drawing.Point(10, 770);
            this.resetButton.Name = "resetButton";
            this.resetButton.Size = new System.Drawing.Size(170, 25);
            this.resetButton.TabIndex = 21;
            this.resetButton.Text = "Reset";
            this.resetButton.FlatStyle = System.Windows.Forms.FlatStyle.Flat;
            this.resetButton.FlatAppearance.BorderSize = 1;
            this.resetButton.FlatAppearance.BorderColor = System.Drawing.Color.Black;
            this.resetButton.BackColor = System.Drawing.SystemColors.Control;
            this.resetButton.ForeColor = System.Drawing.Color.Black;
            this.resetButton.Font = new System.Drawing.Font("Consolas", 9F);
            this.resetButton.Click += new System.EventHandler(this.resetButton_Click);

            // outputConsoleTextBox
            this.outputConsoleTextBox.Location = new System.Drawing.Point(12, 540);
            this.outputConsoleTextBox.Multiline = true;
            this.outputConsoleTextBox.Name = "outputConsoleTextBox";
            this.outputConsoleTextBox.Size = new System.Drawing.Size(860, 150);
            this.outputConsoleTextBox.TabIndex = 22;
            this.outputConsoleTextBox.ReadOnly = true;
            this.outputConsoleTextBox.ScrollBars = System.Windows.Forms.ScrollBars.Vertical;
            this.outputConsoleTextBox.Font = new System.Drawing.Font("Consolas", 9F);
            this.outputConsoleTextBox.BackColor = System.Drawing.Color.White;
            this.outputConsoleTextBox.ForeColor = System.Drawing.Color.Black;
            this.outputConsoleTextBox.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;

            // visualizationPanel
            this.visualizationPanel.Location = new System.Drawing.Point(318, 12);
            this.visualizationPanel.Name = "visualizationPanel";
            this.visualizationPanel.Size = new System.Drawing.Size(560, 520);
            this.visualizationPanel.TabIndex = 23;
            this.visualizationPanel.BackColor = System.Drawing.Color.White;
            this.visualizationPanel.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.visualizationPanel.Paint += new System.Windows.Forms.PaintEventHandler(this.visualizationPanel_Paint);

            // Form1
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(900, 700);
            this.Controls.Add(this.controlPanel);
            this.Controls.Add(this.outputConsoleTextBox);
            this.Controls.Add(this.visualizationPanel);
            this.Name = "Form1";
            this.Text = "Estuarine Circulation Modeling";
            this.Font = new System.Drawing.Font("Consolas", 9F);
            this.BackColor = System.Drawing.SystemColors.Control;
            this.ForeColor = System.Drawing.Color.Black;
            this.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedDialog;
            this.MaximizeBox = false;
            this.controlPanel.ResumeLayout(false);
            this.controlPanel.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();
        }

        private System.Windows.Forms.Panel controlPanel;
        private System.Windows.Forms.Button tidalStrainButton;
        private System.Windows.Forms.Button satButton;
        private System.Windows.Forms.Button compForcingButton;
        private System.Windows.Forms.Button tamButton;
        private System.Windows.Forms.Button atmButton;
        private System.Windows.Forms.Label riverInflowLabel;
        private System.Windows.Forms.TextBox riverInflowTextBox;
        private System.Windows.Forms.Label tidalAmplitudeLabel;
        private System.Windows.Forms.TextBox tidalAmplitudeTextBox;
        private System.Windows.Forms.Label tidalPeriodLabel;
        private System.Windows.Forms.TextBox tidalPeriodTextBox;
        private System.Windows.Forms.Label oceanSalinityLabel;
        private System.Windows.Forms.TextBox oceanSalinityTextBox;
        private System.Windows.Forms.Label oceanTemperatureLabel;
        private System.Windows.Forms.TextBox oceanTemperatureTextBox;
        private System.Windows.Forms.Label estuaryLengthLabel;
        private System.Windows.Forms.TextBox estuaryLengthTextBox;
        private System.Windows.Forms.Label estuaryDepthLabel;
        private System.Windows.Forms.TextBox estuaryDepthTextBox;
        private System.Windows.Forms.CheckBox useRANSSolverCheckBox;
        private System.Windows.Forms.CheckBox useNonHydrostaticCheckBox;
        private System.Windows.Forms.Button startButton;
        private System.Windows.Forms.Button pauseButton;
        private System.Windows.Forms.Button resetButton;
        private System.Windows.Forms.Button stratificationButton;
        private System.Windows.Forms.Button largeEddyButton;
        private System.Windows.Forms.Button lbLesButton;
        private System.Windows.Forms.TextBox outputConsoleTextBox;
        private System.Windows.Forms.Panel visualizationPanel;
    }
}