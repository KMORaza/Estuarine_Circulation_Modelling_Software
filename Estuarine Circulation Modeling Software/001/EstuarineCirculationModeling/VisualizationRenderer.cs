using System;
using System.Drawing;

namespace EstuarineCirculationModeling
{
    public static class VisualizationRenderer
    {
        public static void Render(Graphics g, int width, int height, EstuarineModel model)
        {
            g.Clear(Color.White);

            // Calculate water level
            double waterLevel = model.TidalAmplitude * Math.Sin(2 * Math.PI * model.CurrentTime / model.TidalPeriod);
            int waterLevelPixel = (int)(height / 2 - (waterLevel / model.TidalAmplitude) * (height / 8));

            // Draw estuary background
            using (Brush waterBrush = new SolidBrush(Color.LightBlue))
            {
                g.FillRectangle(waterBrush, 0, waterLevelPixel, width, height - waterLevelPixel);
            }

            // Draw salt wedge
            int wedgePixelX = (int)(model.SaltWedgePosition / model.EstuaryLength * width);
            using (Brush saltBrush = new SolidBrush(Color.DarkBlue))
            {
                g.FillRectangle(saltBrush, wedgePixelX, waterLevelPixel, width - wedgePixelX, height - waterLevelPixel);
            }

            // Draw salinity profile
            using (Pen salinityPen = new Pen(Color.Black, 2))
            {
                Point[] points = new Point[100];
                for (int i = 0; i < 100; i++)
                {
                    double x = i * model.EstuaryLength / 100;
                    double salinity = model.GetSalinityAtPoint(x);
                    int pixelX = i * width / 100;
                    int pixelY = (int)(waterLevelPixel - (salinity / model.SalinityOcean) * (height / 8));
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
                    int pixelX = i * width / 100;
                    int pixelY = (int)(waterLevelPixel - (temperature / model.TemperatureOcean) * (height / 8));
                    points[i] = new Point(pixelX, pixelY);
                }
                g.DrawLines(temperaturePen, points);
            }

            // Draw velocity profile if RANS solver is enabled
            if (model.UseRANSSolver)
            {
                using (Pen velocityPen = new Pen(Color.Red, 1))
                {
                    double[] velocityProfile = model.GetVelocityProfile();
                    int gridPoints = velocityProfile.Length;
                    double maxVelocity = 0.1; // Scale for visualization
                    for (int i = 0; i < gridPoints - 1; i += 5)
                    {
                        int pixelX = i * width / gridPoints;
                        int baseY = height - 20;
                        double velocity = velocityProfile[i];
                        int arrowLength = (int)(velocity / maxVelocity * 30);
                        g.DrawLine(velocityPen, pixelX, baseY, pixelX + arrowLength, baseY);
                        if (arrowLength != 0)
                        {
                            int arrowSize = 5;
                            g.DrawLine(velocityPen, pixelX + arrowLength, baseY, pixelX + arrowLength - arrowSize, baseY - arrowSize);
                            g.DrawLine(velocityPen, pixelX + arrowLength, baseY, pixelX + arrowLength - arrowSize, baseY + arrowSize);
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
                if (model.UseRANSSolver)
                {
                    g.DrawString("Velocity Vectors (Red)", font, Brushes.Red, 10, 50);
                    g.DrawString(model.IsNonHydrostatic() ? "Non-Hydrostatic" : "Hydrostatic", font, Brushes.Purple, 10, 70);
                }
                g.DrawString($"Water Level: {waterLevel:F2}m", font, Brushes.Blue, 10, 90);
            }
        }
    }
}