# Kuang_v0.3.0

## Logging

### v0.4.0

In v0.4.0, I've made following changes:

- Reconstruction of radiative heating is removed from `script/Reconstruct.jl`. Alternatively, it can be diagnosed from the variables predicted by the model and multipled with derived coefficients.

- The reconstruction in `script/Reconstruct.jl` would become a daily-averaged value, so as to decline the usage of disk.
