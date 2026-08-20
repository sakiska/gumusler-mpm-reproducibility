# Pb–Zn–Cu Variogram and Ordinary Kriging Model Comparison

- Lag width: 50 m
- Cutoff distance: 800 m
- Direction: Omnidirectional
- Transformations: Raw and log10
- Candidate models: Spherical, Exponential and Gaussian
- Validation: Leave-one-out Ordinary Kriging cross-validation
- Primary selection criterion: Lowest RMSE on the original ppm scale
- Secondary criteria: Lowest absolute ME and then lowest variogram RSS
- Log10 predictions were back-transformed with an approximate lognormal bias correction.

## Selected combinations

- Pb: Raw + Exponential; RMSE = 165.367 ppm; ME = 0.221 ppm; range = 70.2 m.
- Zn: Raw + Exponential; RMSE = 34.699 ppm; ME = 0.081 ppm; range = 120.8 m.
- Cu: Raw + Gaussian; RMSE = 125.497 ppm; ME = 0 ppm; range = 109.8 m.
