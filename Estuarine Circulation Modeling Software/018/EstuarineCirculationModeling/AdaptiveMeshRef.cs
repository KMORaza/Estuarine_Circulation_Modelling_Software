using System;
using System.Collections.Generic;
using System.Linq; // Added for ToList()

namespace EstuarineCirculationModeling
{
    public class AdaptiveMeshRef
    {
        private readonly int baseGridX = 400;
        private readonly int baseGridZ = 100;
        private readonly double estuaryLength;
        private readonly double estuaryDepth;
        private readonly double[,] salinityProfile;
        private readonly double[,] temperatureProfile;
        private readonly double[,] sedimentProfile;
        private readonly double[,] uVelocityProfile;
        private readonly double[,] wVelocityProfile;
        private readonly double[] bedloadProfile;
        private readonly double kinematicViscosity;
        private readonly double thermalDiffusivity;
        private readonly double sedimentDiffusivity;
        private readonly double settlingVelocity;
        private readonly double coriolisParameter;
        private readonly double turbulenceIntensity;
        private readonly double tidalPeriod;
        private readonly double tidalAmplitude;
        private readonly double windSpeed;
        private readonly double windDirection;
        private readonly double windPeriod;
        private readonly double riverDischargeAmplitude;
        private readonly double riverPeriod;
        private readonly Random rand = new Random();
        private readonly double rho0 = 1000.0;
        private readonly double rhoA = 1.225;
        private readonly double rhoS = 2650.0;
        private readonly double Cd = 0.0012;
        private readonly double alpha = 2e-4;
        private readonly double beta = 7.6e-4;
        private readonly double gamma = (2650.0 - 1000.0) / (1000.0 * 2650.0);
        private readonly double g = 9.81;
        private readonly double d = 0.001;
        private readonly double thetaC = 0.047;
        private readonly double alphaE = 0.0001;
        private readonly double kw = 0.1;
        private readonly double T0 = 15.0;
        private readonly double S0 = 35.0;
        private bool[,] refineFlags; // True for cells needing refinement
        private readonly Dictionary<(int, int), double[,]> fineSalinityGrids;
        private readonly Dictionary<(int, int), double[,]> fineTemperatureGrids;
        private readonly Dictionary<(int, int), double[,]> fineSedimentGrids;
        private readonly Dictionary<(int, int), double[,]> fineUVelocityGrids;
        private readonly Dictionary<(int, int), double[,]> fineWVelocityGrids;
        private readonly int refineFactor = 2; // 2x refinement
        private readonly double gradientThresholdFactor = 0.9; // Refine top 10% of gradients

        public AdaptiveMeshRef(
            double estuaryLength, double estuaryDepth, double tidalPeriod, double tidalAmplitude,
            double[,] salinityProfile, double[,] temperatureProfile, double[,] sedimentProfile,
            double[,] uVelocityProfile, double[,] wVelocityProfile, double[] bedloadProfile,
            double kinematicViscosity, double thermalDiffusivity, double sedimentDiffusivity,
            double settlingVelocity, double coriolisParameter, double turbulenceIntensity,
            double windSpeed, double windDirection, double windPeriod,
            double riverDischargeAmplitude, double riverPeriod)
        {
            this.estuaryLength = estuaryLength;
            this.estuaryDepth = estuaryDepth;
            this.tidalPeriod = tidalPeriod;
            this.tidalAmplitude = tidalAmplitude;
            this.salinityProfile = salinityProfile;
            this.temperatureProfile = temperatureProfile;
            this.sedimentProfile = sedimentProfile;
            this.uVelocityProfile = uVelocityProfile;
            this.wVelocityProfile = wVelocityProfile;
            this.bedloadProfile = bedloadProfile;
            this.kinematicViscosity = kinematicViscosity;
            this.thermalDiffusivity = thermalDiffusivity;
            this.sedimentDiffusivity = sedimentDiffusivity;
            this.settlingVelocity = settlingVelocity;
            this.coriolisParameter = coriolisParameter;
            this.turbulenceIntensity = turbulenceIntensity;
            this.windSpeed = windSpeed;
            this.windDirection = windDirection;
            this.windPeriod = windPeriod;
            this.riverDischargeAmplitude = riverDischargeAmplitude;
            this.riverPeriod = riverPeriod;
            refineFlags = new bool[baseGridX, baseGridZ];
            fineSalinityGrids = new Dictionary<(int, int), double[,]>();
            fineTemperatureGrids = new Dictionary<(int, int), double[,]>();
            fineSedimentGrids = new Dictionary<(int, int), double[,]>();
            fineUVelocityGrids = new Dictionary<(int, int), double[,]>();
            fineWVelocityGrids = new Dictionary<(int, int), double[,]>();
        }

        public void RefineGrid(double currentTime)
        {
            double dx = estuaryLength / (baseGridX - 1);
            double dz = estuaryDepth / (baseGridZ - 1);

            // Compute gradients and find threshold
            List<double> gradientMagnitudes = new List<double>();
            for (int i = 1; i < baseGridX - 1; i++)
            {
                for (int j = 1; j < baseGridZ - 1; j++)
                {
                    // Salinity gradients
                    double dsdx = (salinityProfile[i + 1, j] - salinityProfile[i - 1, j]) / (2 * dx);
                    double dsdz = (salinityProfile[i, j + 1] - salinityProfile[i, j - 1]) / (2 * dz);
                    // Temperature gradients
                    double dTdx = (temperatureProfile[i + 1, j] - temperatureProfile[i - 1, j]) / (2 * dx);
                    double dTdz = (temperatureProfile[i, j + 1] - temperatureProfile[i, j - 1]) / (2 * dz);
                    // Velocity gradients
                    double dudx = (uVelocityProfile[i + 1, j] - uVelocityProfile[i - 1, j]) / (2 * dx);
                    double dudz = (uVelocityProfile[i, j + 1] - uVelocityProfile[i, j - 1]) / (2 * dz);
                    double gradientMag = Math.Sqrt(
                        dsdx * dsdx + dsdz * dsdz +
                        dTdx * dTdx + dTdz * dTdz +
                        dudx * dudx + dudz * dudz);
                    gradientMagnitudes.Add(gradientMag);
                }
            }
            gradientMagnitudes.Sort();
            double threshold = gradientMagnitudes[(int)(gradientMagnitudes.Count * gradientThresholdFactor)];

            // Flag cells for refinement
            refineFlags = new bool[baseGridX, baseGridZ];
            fineSalinityGrids.Clear();
            fineTemperatureGrids.Clear();
            fineSedimentGrids.Clear();
            fineUVelocityGrids.Clear();
            fineWVelocityGrids.Clear();
            for (int i = 1; i < baseGridX - 1; i++)
            {
                for (int j = 1; j < baseGridZ - 1; j++)
                {
                    double dsdx = (salinityProfile[i + 1, j] - salinityProfile[i - 1, j]) / (2 * dx);
                    double dsdz = (salinityProfile[i, j + 1] - salinityProfile[i, j - 1]) / (2 * dz);
                    double dTdx = (temperatureProfile[i + 1, j] - temperatureProfile[i - 1, j]) / (2 * dx);
                    double dTdz = (temperatureProfile[i, j + 1] - temperatureProfile[i, j - 1]) / (2 * dz);
                    double dudx = (uVelocityProfile[i + 1, j] - uVelocityProfile[i - 1, j]) / (2 * dx);
                    double dudz = (uVelocityProfile[i, j + 1] - uVelocityProfile[i, j - 1]) / (2 * dz);
                    double gradientMag = Math.Sqrt(
                        dsdx * dsdx + dsdz * dsdz +
                        dTdx * dTdx + dTdz * dTdz +
                        dudx * dudx + dudz * dudz);
                    if (gradientMag > threshold)
                    {
                        refineFlags[i, j] = true;
                        // Initialize fine grids (2x2 cells per coarse cell)
                        fineSalinityGrids[(i, j)] = new double[refineFactor, refineFactor];
                        fineTemperatureGrids[(i, j)] = new double[refineFactor, refineFactor];
                        fineSedimentGrids[(i, j)] = new double[refineFactor, refineFactor];
                        fineUVelocityGrids[(i, j)] = new double[refineFactor, refineFactor];
                        fineWVelocityGrids[(i, j)] = new double[refineFactor, refineFactor];
                        // Interpolate initial values
                        for (int fi = 0; fi < refineFactor; fi++)
                        {
                            for (int fj = 0; fj < refineFactor; fj++)
                            {
                                double xFrac = fi / (double)refineFactor;
                                double zFrac = fj / (double)refineFactor;
                                fineSalinityGrids[(i, j)][fi, fj] = Interpolate(
                                    salinityProfile[i, j], salinityProfile[i + 1, j],
                                    salinityProfile[i, j + 1], salinityProfile[i + 1, j + 1], xFrac, zFrac);
                                fineTemperatureGrids[(i, j)][fi, fj] = Interpolate(
                                    temperatureProfile[i, j], temperatureProfile[i + 1, j],
                                    temperatureProfile[i, j + 1], temperatureProfile[i + 1, j + 1], xFrac, zFrac);
                                fineSedimentGrids[(i, j)][fi, fj] = Interpolate(
                                    sedimentProfile[i, j], sedimentProfile[i + 1, j],
                                    sedimentProfile[i, j + 1], sedimentProfile[i + 1, j + 1], xFrac, zFrac);
                                fineUVelocityGrids[(i, j)][fi, fj] = Interpolate(
                                    uVelocityProfile[i, j], uVelocityProfile[i + 1, j],
                                    uVelocityProfile[i, j + 1], uVelocityProfile[i + 1, j + 1], xFrac, zFrac);
                                fineWVelocityGrids[(i, j)][fi, fj] = Interpolate(
                                    wVelocityProfile[i, j], wVelocityProfile[i + 1, j],
                                    wVelocityProfile[i, j + 1], wVelocityProfile[i + 1, j + 1], xFrac, zFrac);
                            }
                        }
                    }
                }
            }
        }

        private double Interpolate(double v00, double v10, double v01, double v11, double xFrac, double zFrac)
        {
            // Bilinear interpolation
            return (1 - xFrac) * (1 - zFrac) * v00 +
                   xFrac * (1 - zFrac) * v10 +
                   (1 - xFrac) * zFrac * v01 +
                   xFrac * zFrac * v11;
        }

        public void UpdateFields(double currentTime, double dt)
        {
            double baseDx = estuaryLength / (baseGridX - 1);
            double baseDz = estuaryDepth / (baseGridZ - 1);
            double fineDx = baseDx / refineFactor;
            double fineDz = baseDz / refineFactor;
            double fineDt = dt / (refineFactor * refineFactor); 
            int subSteps = refineFactor * refineFactor;

            // Compute wind and river conditions
            double windPhase = 2 * Math.PI * currentTime / windPeriod;
            double currentWindSpeed = windSpeed * Math.Abs(Math.Sin(windPhase));
            double currentWindDirection = (windDirection + 30.0 * Math.Sin(windPhase)) % 360.0;
            double tauW = rhoA * Cd * currentWindSpeed * currentWindSpeed;
            double tauWx = tauW * Math.Cos(currentWindDirection * Math.PI / 180.0);
            double tauWz = tauW * Math.Sin(currentWindDirection * Math.PI / 180.0);
            double ekmanDepth = Math.Sqrt(2 * kinematicViscosity / Math.Max(coriolisParameter, 1e-6));
            double uEkman = coriolisParameter != 0 ? tauWx / (rho0 * coriolisParameter * ekmanDepth) : 0.0;
            double wEkman = coriolisParameter != 0 ? tauWz / (rho0 * coriolisParameter * ekmanDepth) : 0.0;
            double riverPhase = 2 * Math.PI * currentTime / riverPeriod;
            double baseDischarge = 0.1;
            double currentRiverDischarge = baseDischarge + riverDischargeAmplitude * Math.Sin(riverPhase);
            currentRiverDischarge = Math.Max(0.0, currentRiverDischarge);
            double riverSalinity = 10.0 * (1.0 - currentRiverDischarge / (baseDischarge + riverDischargeAmplitude));
            double riverTemperature = 14.0 - 4.0 * (currentRiverDischarge / (baseDischarge + riverDischargeAmplitude));
            double riverSediment = 200.0 * (currentRiverDischarge / (baseDischarge + riverDischargeAmplitude));
            double tidalPhase = 2 * Math.PI * currentTime / tidalPeriod;
            double tidalVelocity = tidalAmplitude * Math.Cos(tidalPhase);

            // Update refined regions
            var keys = fineSalinityGrids.Keys.ToList(); // Create a snapshot of the keys
            foreach (var key in keys)
            {
                int i = key.Item1;
                int j = key.Item2;
                double[,] newSalinity = new double[refineFactor, refineFactor];
                double[,] newTemperature = new double[refineFactor, refineFactor];
                double[,] newSediment = new double[refineFactor, refineFactor];
                double[,] newUVelocity = new double[refineFactor, refineFactor];
                double[,] newWVelocity = new double[refineFactor, refineFactor];

                // Perform sub-steps
                for (int step = 0; step < subSteps; step++)
                {
                    double[,] tempSalinity = (double[,])fineSalinityGrids[key].Clone();
                    double[,] tempTemperature = (double[,])fineTemperatureGrids[key].Clone();
                    double[,] tempSediment = (double[,])fineSedimentGrids[key].Clone();
                    double[,] tempUVelocity = (double[,])fineUVelocityGrids[key].Clone();
                    double[,] tempWVelocity = (double[,])fineWVelocityGrids[key].Clone();

                    for (int fi = 1; fi < refineFactor - 1; fi++)
                    {
                        for (int fj = 1; fj < refineFactor - 1; fj++)
                        {
                            // Update velocities
                            double dudx = (tempUVelocity[fi + 1, fj] - tempUVelocity[fi - 1, fj]) / (2 * fineDx);
                            double dudz = (tempUVelocity[fi, fj + 1] - tempUVelocity[fi, fj - 1]) / (2 * fineDz);
                            double uAdvection = tempUVelocity[fi, fj] * dudx + tempWVelocity[fi, fj] * dudz;
                            double d2udx2 = (tempUVelocity[fi + 1, fj] - 2 * tempUVelocity[fi, fj] + tempUVelocity[fi - 1, fj]) / (fineDx * fineDx);
                            double d2udz2 = (tempUVelocity[fi, fj + 1] - 2 * tempUVelocity[fi, fj] + tempUVelocity[fi, fj - 1]) / (fineDz * fineDz);
                            double uViscosity = kinematicViscosity * (d2udx2 + d2udz2);
                            double uCoriolis = -coriolisParameter * tempWVelocity[fi, fj];
                            double uForcing = (tidalVelocity - tempUVelocity[fi, fj]) / tidalPeriod;
                            double zFraction = ((double)j + fj / (double)refineFactor) / (baseGridZ - 1);
                            double windMixing = kw * tauW * (1.0 - zFraction);
                            double uWind = fj >= refineFactor - 1 && j >= baseGridZ - 5 ? uEkman * Math.Exp(-zFraction / ekmanDepth) : 0.0;
                            newUVelocity[fi, fj] = tempUVelocity[fi, fj] + fineDt * (-uAdvection + uViscosity + uCoriolis + uForcing + uWind);
                            newUVelocity[fi, fj] += (turbulenceIntensity + windMixing) * tidalAmplitude * (rand.NextDouble() - 0.5);

                            double dwdx = (tempWVelocity[fi + 1, fj] - tempWVelocity[fi - 1, fj]) / (2 * fineDx);
                            double dwdz = (tempWVelocity[fi, fj + 1] - tempWVelocity[fi, fj - 1]) / (2 * fineDz);
                            double wAdvection = tempUVelocity[fi, fj] * dwdx + tempWVelocity[fi, fj] * dwdz;
                            double d2wdx2 = (tempWVelocity[fi + 1, fj] - 2 * tempWVelocity[fi, fj] + tempWVelocity[fi - 1, fj]) / (fineDx * fineDx);
                            double d2wdz2 = (tempWVelocity[fi, fj + 1] - 2 * tempWVelocity[fi, fj] + tempWVelocity[fi, fj - 1]) / (fineDz * fineDz);
                            double wViscosity = kinematicViscosity * (d2wdx2 + d2wdz2);
                            double wCoriolis = coriolisParameter * tempUVelocity[fi, fj];
                            double wWind = fj >= refineFactor - 1 && j >= baseGridZ - 5 ? wEkman * Math.Exp(-zFraction / ekmanDepth) : 0.0;
                            newWVelocity[fi, fj] = tempWVelocity[fi, fj] + fineDt * (-wAdvection + wViscosity + wCoriolis + wWind);
                            newWVelocity[fi, fj] += (turbulenceIntensity + windMixing) * tidalAmplitude * (rand.NextDouble() - 0.5);

                            // Update salinity
                            double dsdx = (tempSalinity[fi + 1, fj] - tempSalinity[fi - 1, fj]) / (2 * fineDx);
                            double dsdz = (tempSalinity[fi, fj + 1] - tempSalinity[fi, fj - 1]) / (2 * fineDz);
                            double sAdvection = tempUVelocity[fi, fj] * dsdx + tempWVelocity[fi, fj] * dsdz;
                            double d2sdx2 = (tempSalinity[fi + 1, fj] - 2 * tempSalinity[fi, fj] + tempSalinity[fi - 1, fj]) / (fineDx * fineDx);
                            double d2sdz2 = (tempSalinity[fi, fj + 1] - 2 * tempSalinity[fi, fj] + tempSalinity[fi, fj - 1]) / (fineDz * fineDz);
                            double sDiffusion = kinematicViscosity * (d2sdx2 + d2sdz2);
                            newSalinity[fi, fj] = tempSalinity[fi, fj] + fineDt * (-sAdvection + sDiffusion);

                            // Update temperature
                            double dTdx = (tempTemperature[fi + 1, fj] - tempTemperature[fi - 1, fj]) / (2 * fineDx);
                            double dTdz = (tempTemperature[fi, fj + 1] - tempTemperature[fi, fj - 1]) / (2 * fineDz);
                            double tAdvection = tempUVelocity[fi, fj] * dTdx + tempWVelocity[fi, fj] * dTdz;
                            double d2Tdx2 = (tempTemperature[fi + 1, fj] - 2 * tempTemperature[fi, fj] + tempTemperature[fi - 1, fj]) / (fineDx * fineDx);
                            double d2Tdz2 = (tempTemperature[fi, fj + 1] - 2 * tempTemperature[fi, fj] + tempTemperature[fi, fj - 1]) / (fineDz * fineDz);
                            double tDiffusion = thermalDiffusivity * (d2Tdx2 + d2Tdz2);
                            newTemperature[fi, fj] = tempTemperature[fi, fj] + fineDt * (-tAdvection + tDiffusion);

                            // Update sediment
                            double dCdx = (tempSediment[fi + 1, fj] - tempSediment[fi - 1, fj]) / (2 * fineDx);
                            double dCdz = (tempSediment[fi, fj + 1] - tempSediment[fi, fj - 1]) / (2 * fineDz);
                            double cAdvection = tempUVelocity[fi, fj] * dCdx + (tempWVelocity[fi, fj] + settlingVelocity) * dCdz;
                            double d2Cdx2 = (tempSediment[fi + 1, fj] - 2 * tempSediment[fi, fj] + tempSediment[fi - 1, fj]) / (fineDx * fineDx);
                            double d2Cdz2 = (tempSediment[fi, fj + 1] - 2 * tempSediment[fi, fj] + tempSediment[fi, fj - 1]) / (fineDz * fineDz);
                            double cDiffusion = sedimentDiffusivity * (d2Cdx2 + d2Cdz2);
                            newSediment[fi, fj] = Math.Max(0.0, tempSediment[fi, fj] + fineDt * (-cAdvection + cDiffusion));
                        }
                    }

                    // Apply boundary conditions within fine grid
                    for (int fi = 0; fi < refineFactor; fi++)
                    {
                        for (int fj = 0; fj < refineFactor; fj++)
                        {
                            bool isBoundary = fi == 0 || fi == refineFactor - 1 || fj == 0 || fj == refineFactor - 1;
                            if (isBoundary)
                            {
                                double xFrac = fi / (double)refineFactor;
                                double zFrac = fj / (double)refineFactor;
                                newSalinity[fi, fj] = Interpolate(
                                    salinityProfile[i, j], salinityProfile[i + 1, j],
                                    salinityProfile[i, j + 1], salinityProfile[i + 1, j + 1], xFrac, zFrac);
                                newTemperature[fi, fj] = Interpolate(
                                    temperatureProfile[i, j], temperatureProfile[i + 1, j],
                                    temperatureProfile[i, j + 1], temperatureProfile[i + 1, j + 1], xFrac, zFrac);
                                newSediment[fi, fj] = Interpolate(
                                    sedimentProfile[i, j], sedimentProfile[i + 1, j],
                                    sedimentProfile[i, j + 1], sedimentProfile[i + 1, j + 1], xFrac, zFrac);
                                newUVelocity[fi, fj] = Interpolate(
                                    uVelocityProfile[i, j], uVelocityProfile[i + 1, j],
                                    uVelocityProfile[i, j + 1], uVelocityProfile[i + 1, j + 1], xFrac, zFrac);
                                newWVelocity[fi, fj] = Interpolate(
                                    wVelocityProfile[i, j], wVelocityProfile[i + 1, j],
                                    wVelocityProfile[i, j + 1], wVelocityProfile[i + 1, j + 1], xFrac, zFrac);
                            }
                            if (i == 0) // River boundary
                            {
                                newUVelocity[0, fj] = currentRiverDischarge;
                                newWVelocity[0, fj] = 0.0;
                                newSalinity[0, fj] = riverSalinity;
                                newTemperature[0, fj] = riverTemperature;
                                newSediment[0, fj] = riverSediment;
                            }
                            if (i == baseGridX - 2) // Ocean boundary
                            {
                                newUVelocity[refineFactor - 1, fj] = tidalVelocity;
                                newWVelocity[refineFactor - 1, fj] = 0.0;
                                newSalinity[refineFactor - 1, fj] = 35.0;
                                newTemperature[refineFactor - 1, fj] = 15.0;
                                newSediment[refineFactor - 1, fj] = 0.0;
                            }
                            if (j == 0) // Bottom
                            {
                                newUVelocity[fi, 0] = 0.0;
                                newWVelocity[fi, 0] = 0.0;
                                newSalinity[fi, 0] = newSalinity[fi, 1];
                                newTemperature[fi, 0] = newTemperature[fi, 1];
                                newSediment[fi, 0] = newSediment[fi, 1];
                            }
                            if (j == baseGridZ - 2) // Surface
                            {
                                newUVelocity[fi, refineFactor - 1] = tidalVelocity + tauWx / (rho0 * kinematicViscosity);
                                newWVelocity[fi, refineFactor - 1] = tauWz / (rho0 * kinematicViscosity);
                                newSalinity[fi, refineFactor - 1] = newSalinity[fi, refineFactor - 2];
                                newTemperature[fi, refineFactor - 1] = 15.0;
                                newSediment[fi, refineFactor - 1] = newSediment[fi, refineFactor - 2];
                            }
                        }
                    }
                    fineSalinityGrids[key] = (double[,])newSalinity.Clone();
                    fineTemperatureGrids[key] = (double[,])newTemperature.Clone();
                    fineSedimentGrids[key] = (double[,])newSediment.Clone();
                    fineUVelocityGrids[key] = (double[,])newUVelocity.Clone();
                    fineWVelocityGrids[key] = (double[,])newWVelocity.Clone();
                }
            }

            // Update coarse grid with fine grid values (average)
            foreach (var key in fineSalinityGrids.Keys)
            {
                int i = key.Item1;
                int j = key.Item2;
                double sSum = 0.0, tSum = 0.0, cSum = 0.0, uSum = 0.0, wSum = 0.0;
                for (int fi = 0; fi < refineFactor; fi++)
                {
                    for (int fj = 0; fj < refineFactor; fj++)
                    {
                        sSum += fineSalinityGrids[key][fi, fj];
                        tSum += fineTemperatureGrids[key][fi, fj];
                        cSum += fineSedimentGrids[key][fi, fj];
                        uSum += fineUVelocityGrids[key][fi, fj];
                        wSum += fineWVelocityGrids[key][fi, fj];
                    }
                }
                int count = refineFactor * refineFactor;
                salinityProfile[i, j] = sSum / count;
                temperatureProfile[i, j] = tSum / count;
                sedimentProfile[i, j] = cSum / count;
                uVelocityProfile[i, j] = uSum / count;
                wVelocityProfile[i, j] = wSum / count;
            }

            // Update bedload on coarse grid
            double[] erosionRate = new double[baseGridX];
            double[] depositionRate = new double[baseGridX];
            for (int i = 1; i < baseGridX - 1; i++)
            {
                double dudz = (uVelocityProfile[i, 1] - uVelocityProfile[i, 0]) / baseDz;
                double tauB = rho0 * kinematicViscosity * dudz;
                double uStar = Math.Sqrt(Math.Abs(tauB) / rho0);
                erosionRate[i] = alphaE * uStar * uStar;
                depositionRate[i] = -settlingVelocity * sedimentProfile[i, 1];
            }
            for (int i = 1; i < baseGridX - 1; i++)
            {
                double dudz = (uVelocityProfile[i, 1] - uVelocityProfile[i, 0]) / baseDz;
                double tauB = rho0 * kinematicViscosity * dudz;
                double theta = tauB / ((rhoS - rho0) * g * d);
                double qb = theta > thetaC ? 8.0 * Math.Pow(theta - thetaC, 1.5) * Math.Sqrt((rhoS / rho0 - 1) * g * d * d * d) : 0.0;
                double qbNext = 0.0;
                if (i < baseGridX - 2)
                {
                    double dudzNext = (uVelocityProfile[i + 1, 1] - uVelocityProfile[i + 1, 0]) / baseDz;
                    double tauBNext = rho0 * kinematicViscosity * dudzNext;
                    double thetaNext = tauBNext / ((rhoS - rho0) * g * d);
                    qbNext = thetaNext > thetaC ? 8.0 * Math.Pow(thetaNext - thetaC, 1.5) * Math.Sqrt((rhoS / rho0 - 1) * g * d * d * d) : 0.0;
                }
                double dqb_dx = (qbNext - qb) / baseDx;
                bedloadProfile[i] = Math.Max(0.0, bedloadProfile[i] + dt * (erosionRate[i] - depositionRate[i] - dqb_dx));
            }
        }

        public bool IsRefined(int i, int j)
        {
            return refineFlags[i, j];
        }

        public double GetSalinity(int i, int j, int fi, int fj)
        {
            return fineSalinityGrids.ContainsKey((i, j)) ? fineSalinityGrids[(i, j)][fi, fj] : salinityProfile[i, j];
        }

        public double GetTemperature(int i, int j, int fi, int fj)
        {
            return fineTemperatureGrids.ContainsKey((i, j)) ? fineTemperatureGrids[(i, j)][fi, fj] : temperatureProfile[i, j];
        }

        public double GetSediment(int i, int j, int fi, int fj)
        {
            return fineSedimentGrids.ContainsKey((i, j)) ? fineSedimentGrids[(i, j)][fi, fj] : sedimentProfile[i, j];
        }

        public double GetUVelocity(int i, int j, int fi, int fj)
        {
            return fineUVelocityGrids.ContainsKey((i, j)) ? fineUVelocityGrids[(i, j)][fi, fj] : uVelocityProfile[i, j];
        }

        public double GetWVelocity(int i, int j, int fi, int fj)
        {
            return fineWVelocityGrids.ContainsKey((i, j)) ? fineWVelocityGrids[(i, j)][fi, fj] : wVelocityProfile[i, j];
        }
    }
}