using System;
using System.Drawing;

namespace EstuarineCirculationModeling
{
    public static class VisualizationRenderer
    {
        public static void Render(Graphics g, int width, int height, EstuarineModel model, Stratification stratification, PassiveScalarTransportEq passiveScalarTransport, double riverScalarConcentration)
        {
            g.Clear(Color.White);

            // Calculate water level
            double waterLevel = model.TidalAmplitude * Math.Sin(2 * Math.PI * model.CurrentTime / model.TidalPeriod);
            int waterLevelPixel = (int)(height / 2 - (waterLevel / Math.Max(0.1, model.TidalAmplitude)) * (height / 8));
            waterLevelPixel = Math.Max(0, Math.Min(height, waterLevelPixel)); // Clamp to panel bounds

            // Draw estuary background
            using (Brush waterBrush = new SolidBrush(Color.LightBlue))
            {
                g.FillRectangle(waterBrush, 0, waterLevelPixel, width, Math.Max(0, height - waterLevelPixel));
            }

            // Draw salt wedge
            int wedgePixelX = (int)(model.SaltWedgePosition / model.EstuaryLength * width);
            wedgePixelX = Math.Max(0, Math.Min(width, wedgePixelX)); // Clamp to panel bounds
            using (Brush saltBrush = new SolidBrush(Color.DarkBlue))
            {
                g.FillRectangle(saltBrush, wedgePixelX, waterLevelPixel, Math.Max(0, width - wedgePixelX), Math.Max(0, height - waterLevelPixel));
            }

            // Draw salinity profile
            using (Pen salinityPen = new Pen(Color.Black, 2))
            {
                Point[] points = new Point[100];
                for (int i = 0; i < 100; i++)
                {
                    double x = i * model.EstuaryLength / 100;
                    double salinity = model.GetSalinityAtPoint(x);
                    salinity = Math.Max(0.0, Math.Min(model.SalinityOcean, salinity)); // Clamp salinity
                    int pixelX = i * width / 100;
                    pixelX = Math.Max(0, Math.Min(width, pixelX)); // Clamp x
                    double salinityRatio = model.SalinityOcean > 0.1 ? salinity / model.SalinityOcean : 0.0; // Prevent division by zero
                    int pixelY = (int)(waterLevelPixel - salinityRatio * (height / 8));
                    pixelY = Math.Max(0, Math.Min(height, pixelY)); // Clamp y
                    points[i] = new Point(pixelX, pixelY);
                }
                g.DrawLines(salinityPen, points);
            }

            // Draw temperature profile
            using (Pen temperaturePen = new Pen(Color.Orange, 2))
            {
                Point[] points = new Point[100];
                for (int i = 0; i < 100; i++)
                {
                    double x = i * model.EstuaryLength / 100;
                    double temperature = model.GetTemperatureAtPoint(x);
                    temperature = Math.Max(0.0, Math.Min(model.TemperatureOcean, temperature)); // Clamp temperature
                    int pixelX = i * width / 100;
                    pixelX = Math.Max(0, Math.Min(width, pixelX)); // Clamp x
                    double temperatureRatio = model.TemperatureOcean > 0.1 ? temperature / model.TemperatureOcean : 0.0; // Prevent division by zero
                    int pixelY = (int)(waterLevelPixel - temperatureRatio * (height / 8));
                    pixelY = Math.Max(0, Math.Min(height, pixelY)); // Clamp y
                    points[i] = new Point(pixelX, pixelY);
                }
                g.DrawLines(temperaturePen, points);
            }

            // Draw passive scalar profile
            using (Pen scalarPen = new Pen(Color.Red, 2))
            {
                double[] passiveScalar = stratification.GetPassiveScalar();
                Point[] points = new Point[100];
                double maxScalar = Math.Max(1.0, riverScalarConcentration); // Adjust for user input
                for (int i = 0; i < 100; i++)
                {
                    double x = i * model.EstuaryLength / 100;
                    int index = (int)(x / model.EstuaryLength * 100);
                    index = Math.Max(0, Math.Min(99, index));
                    double scalar = passiveScalar[index];
                    scalar = Math.Max(0.0, Math.Min(maxScalar, scalar)); // Clamp scalar
                    int pixelX = i * width / 100;
                    pixelX = Math.Max(0, Math.Min(width, pixelX)); // Clamp x
                    double scalarRatio = maxScalar > 0.1 ? scalar / maxScalar : 0.0; // Prevent division by zero
                    int pixelY = (int)(waterLevelPixel - scalarRatio * (height / 8));
                    pixelY = Math.Max(0, Math.Min(height, pixelY)); // Clamp y
                    points[i] = new Point(pixelX, pixelY);
                }
                g.DrawLines(scalarPen, points);
            }

            // Draw velocity profile if RANS solver is enabled
            if (model.UseRANSSolver)
            {
                using (Pen velocityPen = new Pen(Color.Green, 1)) // Changed to Green to avoid confusion with scalar
                {
                    double[] velocityProfile = model.GetVelocityProfile();
                    int gridPoints = velocityProfile.Length;
                    double maxVelocity = 0.1; // Scale for visualization
                    for (int i = 0; i < gridPoints - 1; i += 5)
                    {
                        int pixelX = i * width / gridPoints;
                        pixelX = Math.Max(0, Math.Min(width, pixelX)); // Clamp x
                        int baseY = height - 20;
                        baseY = Math.Max(0, Math.Min(height, baseY)); // Clamp y
                        double velocity = Math.Max(-maxVelocity, Math.Min(maxVelocity, velocityProfile[i])); // Clamp velocity
                        int arrowLength = (int)(velocity / maxVelocity * 30);
                        int endX = pixelX + arrowLength;
                        endX = Math.Max(0, Math.Min(width, endX)); // Clamp endX
                        g.DrawLine(velocityPen, pixelX, baseY, endX, baseY);
                        if (arrowLength != 0)
                        {
                            int arrowSize = 5;
                            int arrowTipY1 = baseY - arrowSize;
                            int arrowTipY2 = baseY + arrowSize;
                            arrowTipY1 = Math.Max(0, Math.Min(height, arrowTipY1)); // Clamp arrow tips
                            arrowTipY2 = Math.Max(0, Math.Min(height, arrowTipY2));
                            g.DrawLine(velocityPen, endX, baseY, endX - arrowSize, arrowTipY1);
                            g.DrawLine(velocityPen, endX, baseY, endX - arrowSize, arrowTipY2);
                        }
                    }
                }
            }

            // Draw labels using Tahoma font
            using (Font font = new Font("Tahoma", 8.25F))
            {
                g.DrawString("River", font, Brushes.Black, 10, height / 2);
                g.DrawString("Ocean", font, Brushes.Black, width - 50, height / 2);
                g.DrawString("Salinity Profile (Black)", font, Brushes.Black, 10, 10);
                g.DrawString("Temperature Profile (Orange)", font, Brushes.Orange, 10, 30);
                g.DrawString("Passive Scalar Profile (Red)", font, Brushes.Red, 10, 50);
                if (model.UseRANSSolver)
                {
                    g.DrawString("Velocity Vectors (Green)", font, Brushes.Green, 10, 70);
                    g.DrawString(model.IsNonHydrostatic() ? "Non-Hydrostatic" : "Hydrostatic", font, Brushes.Purple, 10, 90);
                }
                g.DrawString($"Water Level: {waterLevel:F2}m", font, Brushes.Blue, 10, 110);
            }
        }
    }
}