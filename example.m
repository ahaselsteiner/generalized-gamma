rng('default'); % For reproduceability.

pdTrue = GeneralizedGamma(1.4, 1.0, 1.6);
n = 10000;
sample = pdTrue.drawSample(n);

pdEstimated = GeneralizedGamma();
pdEstimated.fitDist(sample)
