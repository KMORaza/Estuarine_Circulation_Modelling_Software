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
            this.controlPanel = new System.Windows.Forms.GroupBox();
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
            this.useRANSSolverCheckBox = new System.Windows.Forms.CheckBox();
            this.useNonHydrostaticCheckBox = new System.Windows.Forms.CheckBox();
            this.startButton = new System.Windows.Forms.Button();
            this.pauseButton = new System.Windows.Forms.Button();
            this.resetButton = new System.Windows.Forms.Button();
            this.outputConsoleTextBox = new System.Windows.Forms.TextBox();
            this.visualizationPanel = new System.Windows.Forms.Panel();
            this.controlPanel.SuspendLayout();
            this.SuspendLayout();

            // controlPanel
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
            this.controlPanel.Size = new System.Drawing.Size(250, 450); // Increased height for new controls
            this.controlPanel.TabIndex = 0;
            this.controlPanel.TabStop = false;
            this.controlPanel.Text = "Control Panel";
            this.controlPanel.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // riverInflowLabel
            this.riverInflowLabel.AutoSize = true;
            this.riverInflowLabel.Location = new System.Drawing.Point(10, 20);
            this.riverInflowLabel.Name = "riverInflowLabel";
            this.riverInflowLabel.Size = new System.Drawing.Size(80, 13);
            this.riverInflowLabel.TabIndex = 0;
            this.riverInflowLabel.Text = "River Inflow (m³/s):";
            this.riverInflowLabel.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // riverInflowTextBox
            this.riverInflowTextBox.Location = new System.Drawing.Point(10, 40);
            this.riverInflowTextBox.Name = "riverInflowTextBox";
            this.riverInflowTextBox.Size = new System.Drawing.Size(200, 20);
            this.riverInflowTextBox.TabIndex = 1;
            this.riverInflowTextBox.Text = "0.1";
            this.riverInflowTextBox.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // tidalAmplitudeLabel
            this.tidalAmplitudeLabel.AutoSize = true;
            this.tidalAmplitudeLabel.Location = new System.Drawing.Point(10, 70);
            this.tidalAmplitudeLabel.Name = "tidalAmplitudeLabel";
            this.tidalAmplitudeLabel.Size = new System.Drawing.Size(80, 13);
            this.tidalAmplitudeLabel.TabIndex = 2;
            this.tidalAmplitudeLabel.Text = "Tidal Amplitude (m):";
            this.tidalAmplitudeLabel.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // tidalAmplitudeTextBox
            this.tidalAmplitudeTextBox.Location = new System.Drawing.Point(10, 90);
            this.tidalAmplitudeTextBox.Name = "tidalAmplitudeTextBox";
            this.tidalAmplitudeTextBox.Size = new System.Drawing.Size(200, 20);
            this.tidalAmplitudeTextBox.TabIndex = 3;
            this.tidalAmplitudeTextBox.Text = "1.0";
            this.tidalAmplitudeTextBox.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // tidalPeriodLabel
            this.tidalPeriodLabel.AutoSize = true;
            this.tidalPeriodLabel.Location = new System.Drawing.Point(10, 120);
            this.tidalPeriodLabel.Name = "tidalPeriodLabel";
            this.tidalPeriodLabel.Size = new System.Drawing.Size(80, 13);
            this.tidalPeriodLabel.TabIndex = 4;
            this.tidalPeriodLabel.Text = "Tidal Period (s):";
            this.tidalPeriodLabel.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // tidalPeriodTextBox
            this.tidalPeriodTextBox.Location = new System.Drawing.Point(10, 140);
            this.tidalPeriodTextBox.Name = "tidalPeriodTextBox";
            this.tidalPeriodTextBox.Size = new System.Drawing.Size(200, 20);
            this.tidalPeriodTextBox.TabIndex = 5;
            this.tidalPeriodTextBox.Text = "43200";
            this.tidalPeriodTextBox.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // oceanSalinityLabel
            this.oceanSalinityLabel.AutoSize = true;
            this.oceanSalinityLabel.Location = new System.Drawing.Point(10, 170);
            this.oceanSalinityLabel.Name = "oceanSalinityLabel";
            this.oceanSalinityLabel.Size = new System.Drawing.Size(80, 13);
            this.oceanSalinityLabel.TabIndex = 6;
            this.oceanSalinityLabel.Text = "Ocean Salinity (PSU):";
            this.oceanSalinityLabel.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // oceanSalinityTextBox
            this.oceanSalinityTextBox.Location = new System.Drawing.Point(10, 190);
            this.oceanSalinityTextBox.Name = "oceanSalinityTextBox";
            this.oceanSalinityTextBox.Size = new System.Drawing.Size(200, 20);
            this.oceanSalinityTextBox.TabIndex = 7;
            this.oceanSalinityTextBox.Text = "35.0";
            this.oceanSalinityTextBox.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // oceanTemperatureLabel
            this.oceanTemperatureLabel.AutoSize = true;
            this.oceanTemperatureLabel.Location = new System.Drawing.Point(10, 220);
            this.oceanTemperatureLabel.Name = "oceanTemperatureLabel";
            this.oceanTemperatureLabel.Size = new System.Drawing.Size(100, 13);
            this.oceanTemperatureLabel.TabIndex = 8;
            this.oceanTemperatureLabel.Text = "Ocean Temperature (°C):";
            this.oceanTemperatureLabel.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // oceanTemperatureTextBox
            this.oceanTemperatureTextBox.Location = new System.Drawing.Point(10, 240);
            this.oceanTemperatureTextBox.Name = "oceanTemperatureTextBox";
            this.oceanTemperatureTextBox.Size = new System.Drawing.Size(200, 20);
            this.oceanTemperatureTextBox.TabIndex = 9;
            this.oceanTemperatureTextBox.Text = "20.0";
            this.oceanTemperatureTextBox.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // estuaryLengthLabel
            this.estuaryLengthLabel.AutoSize = true;
            this.estuaryLengthLabel.Location = new System.Drawing.Point(10, 270);
            this.estuaryLengthLabel.Name = "estuaryLengthLabel";
            this.estuaryLengthLabel.Size = new System.Drawing.Size(80, 13);
            this.estuaryLengthLabel.TabIndex = 10;
            this.estuaryLengthLabel.Text = "Estuary Length (m):";
            this.estuaryLengthLabel.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // estuaryLengthTextBox
            this.estuaryLengthTextBox.Location = new System.Drawing.Point(10, 290);
            this.estuaryLengthTextBox.Name = "estuaryLengthTextBox";
            this.estuaryLengthTextBox.Size = new System.Drawing.Size(200, 20);
            this.estuaryLengthTextBox.TabIndex = 11;
            this.estuaryLengthTextBox.Text = "10000";
            this.estuaryLengthTextBox.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // estuaryDepthLabel
            this.estuaryDepthLabel.AutoSize = true;
            this.estuaryDepthLabel.Location = new System.Drawing.Point(10, 320);
            this.estuaryDepthLabel.Name = "estuaryDepthLabel";
            this.estuaryDepthLabel.Size = new System.Drawing.Size(80, 13);
            this.estuaryDepthLabel.TabIndex = 12;
            this.estuaryDepthLabel.Text = "Estuary Depth (m):";
            this.estuaryDepthLabel.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // estuaryDepthTextBox
            this.estuaryDepthTextBox.Location = new System.Drawing.Point(10, 340);
            this.estuaryDepthTextBox.Name = "estuaryDepthTextBox";
            this.estuaryDepthTextBox.Size = new System.Drawing.Size(200, 20);
            this.estuaryDepthTextBox.TabIndex = 13;
            this.estuaryDepthTextBox.Text = "10.0";
            this.estuaryDepthTextBox.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // useRANSSolverCheckBox
            this.useRANSSolverCheckBox.AutoSize = true;
            this.useRANSSolverCheckBox.Location = new System.Drawing.Point(10, 370);
            this.useRANSSolverCheckBox.Name = "useRANSSolverCheckBox";
            this.useRANSSolverCheckBox.Size = new System.Drawing.Size(120, 17);
            this.useRANSSolverCheckBox.TabIndex = 14;
            this.useRANSSolverCheckBox.Text = "Use RANS Solver";
            this.useRANSSolverCheckBox.Font = new System.Drawing.Font("Tahoma", 8.25F);
            this.useRANSSolverCheckBox.CheckedChanged += new System.EventHandler(this.useRANSSolverCheckBox_CheckedChanged);

            // useNonHydrostaticCheckBox
            this.useNonHydrostaticCheckBox.AutoSize = true;
            this.useNonHydrostaticCheckBox.Location = new System.Drawing.Point(10, 390);
            this.useNonHydrostaticCheckBox.Name = "useNonHydrostaticCheckBox";
            this.useNonHydrostaticCheckBox.Size = new System.Drawing.Size(140, 17);
            this.useNonHydrostaticCheckBox.TabIndex = 15;
            this.useNonHydrostaticCheckBox.Text = "Use Non-Hydrostatic";
            this.useNonHydrostaticCheckBox.Enabled = false;
            this.useNonHydrostaticCheckBox.Font = new System.Drawing.Font("Tahoma", 8.25F);
            this.useNonHydrostaticCheckBox.CheckedChanged += new System.EventHandler(this.useNonHydrostaticCheckBox_CheckedChanged);

            // startButton
            this.startButton.Location = new System.Drawing.Point(10, 410);
            this.startButton.Name = "startButton";
            this.startButton.Size = new System.Drawing.Size(75, 23);
            this.startButton.TabIndex = 16;
            this.startButton.Text = "Start";
            this.startButton.UseVisualStyleBackColor = true;
            this.startButton.Font = new System.Drawing.Font("Tahoma", 8.25F);
            this.startButton.Click += new System.EventHandler(this.startButton_Click);

            // pauseButton
            this.pauseButton.Location = new System.Drawing.Point(90, 410);
            this.pauseButton.Name = "pauseButton";
            this.pauseButton.Size = new System.Drawing.Size(75, 23);
            this.pauseButton.TabIndex = 17;
            this.pauseButton.Text = "Pause";
            this.pauseButton.UseVisualStyleBackColor = true;
            this.pauseButton.Font = new System.Drawing.Font("Tahoma", 8.25F);
            this.pauseButton.Click += new System.EventHandler(this.pauseButton_Click);

            // resetButton
            this.resetButton.Location = new System.Drawing.Point(170, 410);
            this.resetButton.Name = "resetButton";
            this.resetButton.Size = new System.Drawing.Size(75, 23);
            this.resetButton.TabIndex = 18;
            this.resetButton.Text = "Reset";
            this.resetButton.UseVisualStyleBackColor = true;
            this.resetButton.Font = new System.Drawing.Font("Tahoma", 8.25F);
            this.resetButton.Click += new System.EventHandler(this.resetButton_Click);

            // outputConsoleTextBox
            this.outputConsoleTextBox.Location = new System.Drawing.Point(12, 468);
            this.outputConsoleTextBox.Multiline = true;
            this.outputConsoleTextBox.Name = "outputConsoleTextBox";
            this.outputConsoleTextBox.Size = new System.Drawing.Size(776, 150);
            this.outputConsoleTextBox.TabIndex = 19;
            this.outputConsoleTextBox.ReadOnly = true;
            this.outputConsoleTextBox.ScrollBars = ScrollBars.Vertical;
            this.outputConsoleTextBox.Font = new System.Drawing.Font("Tahoma", 8.25F);

            // visualizationPanel
            this.visualizationPanel.Location = new System.Drawing.Point(278, 12);
            this.visualizationPanel.Name = "visualizationPanel";
            this.visualizationPanel.Size = new System.Drawing.Size(510, 450); // Increased height to match controlPanel
            this.visualizationPanel.TabIndex = 20;
            this.visualizationPanel.Paint += new System.Windows.Forms.PaintEventHandler(this.visualizationPanel_Paint);

            // Form1
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(800, 630); // Increased height for larger controlPanel
            this.Controls.Add(this.controlPanel);
            this.Controls.Add(this.outputConsoleTextBox);
            this.Controls.Add(this.visualizationPanel);
            this.Name = "Form1";
            this.Text = "Estuarine Circulation Modeling";
            this.Font = new System.Drawing.Font("Tahoma", 8.25F);
            this.controlPanel.ResumeLayout(false);
            this.controlPanel.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();
        }

        private System.Windows.Forms.GroupBox controlPanel;
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
        private System.Windows.Forms.TextBox outputConsoleTextBox;
        private System.Windows.Forms.Panel visualizationPanel;
    }
}