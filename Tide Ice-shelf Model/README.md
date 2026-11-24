# Rutford tidal-response model

A standalone MATLAB example that simulates Rutford displacement under **M2, S2, K1 and O1** forcing and compares the results with GPS observations. The default entry point reproduces the four-constituent figure from the original Rutford program, displaying the full period from **2005-03-01 to 2005-05-01**. The model is recalculated when the example is run, using observations that have already been prepared. Data-processing scripts and `s_tide` are not included.

## Run the example

Open **`RUN.m`** in MATLAB and click Run. No internet connection, full CATS2008 grid or additional toolboxes are required. The original README records successful execution in MATLAB R2024b. The entire folder can be moved; the entry script locates its files automatically.

Outputs are written to `output/`. Rerunning the example overwrites files with the same names:

- `Rutford_comparison.png`, `.pdf`, `.fig`: a raster image, a vector figure and an editable MATLAB figure.
- `Rutford_comparison.csv`: observation times, observed displacement `observed_m` and simulated displacement after the time shift and interpolation to those times, `simulated_shifted_m`, in m.
- `Rutford_simulation.csv`: the complete simulation, including unshifted time `model_time`, display time `display_time` and displacement `simulated_m`.
- `Rutford_result.mat`: forward-model results, the comparison table, the complete simulation, parameters and tidal harmonic constants.

## Settings used to reproduce the supplied original figure

The default entry point follows the four-constituent section of `ManualFit/Rutford/load_data_scale_comparison.m` and the plotting settings in `plot_scale_comparison.m` from the original project. It uses the original `Drygalski_forward_model.m`, retaining the **additional 300 days of integration**, linear detrending over the entire simulation period and removal of the mean within the display window. The total simulation duration is approximately **360.99 days**, so the shifted series still covers the full March–April period.

The original figure uses the following time mapping:

```matlab
settings.simulation_shift_days = -4.3;
time_sim = datetime(2005,3,1) + days(result.time_days - 4.3);
```

A negative value shifts the series to the left. This example uses the original program's **-4.3-day** shift rather than the -12-day shift used in the earlier revised version. The two models differ in their phase conventions and calculations, so the same shift cannot be applied directly to both. The shift can be changed at the top of `RUN.m`. It represents manual visual alignment and is recorded in the output parameters and the FIG file's UserData. The legend retains the concise wording of the original figure; the shift is not interpreted as a physical response lag.

The observations are the original program's saved `up3 - St` series, packaged from `scale_comparison_results.mat` into `data/Rutford_original_observations.mat`, in cm. This series already includes the background removal performed by the earlier program. **Running this package does not call `s_tide`.** To reproduce the original figure, the default entry point retains the original program's way of connecting observation points, without reprocessing missing observations or the observational background.

Historical verification recorded in the original README (2026-09-23): after restoring MATLAB's default search path and running the example independently, the maximum difference from the complete original `short_cm` series was approximately **2.25 × 10^-10 cm**, and the maximum difference in simulation time was approximately **3.81 × 10^-9 s**. The observation curves matched, and the simulation covered every observation time.

## CATS2008 harmonic constants and reference point

The official name is **CATS2008** (Circum-Antarctic Tidal Simulation version 2008), not CAST2008. Ocean-tide harmonic amplitudes and Greenwich phase lags come from the project's existing extraction from the original `grid_CATS2008` / `hf.CATS2008.out`, not CATS2008_v2023. Only the four constituent constants are included; neither the model grid nor the extraction workflow is packaged here.

Source: Howard, S. L., Erofeeva, S., & Padman, L. (2019), *CATS2008: Circum-Antarctic Tidal Simulation version 2008*, USAP Data Center, [DOI: 10.15784/601235](https://doi.org/10.15784/601235).

| Purpose of location | Latitude (positive north) | Longitude (positive east) |
| --- | ---: | ---: |
| **Ocean-tide forcing reference point** | **-78.800000°** | **-82.000000°** |
| Rutford SEI1 GPS, approximate initial position in December 2004 | -78.141646342° | -88.912341810° |

The GPS station is on grounded ice and is not a valid ocean extraction point in the original CATS2008 model. The example therefore retains the separate ocean-tide reference point selected for the project. This reference point is not the GPS station location, and its tide heights do not represent local tide heights at the GPS station.

| Constituent | Period (h) | Amplitude (m) | Greenwich phase lag (°) |
| --- | ---: | ---: | ---: |
| M2 | 12.4206014422 | 1.5963720089 | 290.3400421509 |
| S2 | 12.0000003575 | 1.0871556105 | 259.1790103835 |
| K1 | 23.9344658347 | 0.4123592524 | 291.6604258229 |
| O1 | 25.8193432501 | 0.4154857471 | 301.6375998135 |

Full precision is retained in `data/Rutford_CATS2008.csv`. Periods and angular frequencies follow the corresponding constituent definitions in TMD. Only these four constituents are supplied as inputs; MSF, MM and minor constituents are not added. The nonlinear model can itself generate long-period responses.

**Phase convention used to reproduce the original figure:** the original script directly uses the numerical CATS2008 phase lags as initial phases of sine functions, then aligns the simulation by shifting the display time. `RUN.m` retains this convention to reproduce the original figure. This is not described as an absolute-time tide-height prediction with astronomical corrections. A Greenwich phase lag is not the same as the initial phase of a sine function at the start of the observation period.

## Model parameters

| Parameter | Value |
| --- | ---: |
| alpha | 1.54 |
| sigma_b0 | 75,000 Pa |
| sigma_h0 | 55,000 Pa |
| sigma_other | 160,000 Pa |
| Length L | 20,000 m |
| Elastic modulus E | 9 × 10^9 Pa |
| Glen parameter A, exponent n=3 | 7 × 10^-25 Pa^-3 s^-1 |
| tidal_amp, multiplier for normalized forcing | 1 |
| Integration time step | 300 s |

The physical parameters are retained from the earlier Rutford case; **they are not supplied by CATS2008**. CATS2008 supplies ocean-tide amplitudes and phases. The model divides the amplitudes by their sum, preserving their relative weights. The original script converts amplitudes to cm before passing them to the model; this unit conversion cancels during normalization. No automatic parameter fitting was performed. Observation values do not enter the model calculations: only the number of observation rows is used to determine the base simulation duration, to which the original core adds 300 days.

`core/Drygalski_forward_model.m` is a complete copy of the original file, retaining its stress-derivative approximation, preliminary integration and numerical integration method. The 5th element of the original version's parameter vector is length L; in the revised version, the 5th element is a delay in hours. The two parameter vectors are not interchangeable.

## Observations and third-party sources

Rutford data: Smith, A., Murray, T., & King, M. (2020), *Long-duration GPS record from Rutford Ice Stream, West Antarctica, 2004–2007 (Version 1.0)*, UK Polar Data Centre, [DOI: 10.5285/dac20505-a56e-4beb-97ba-077eecd587c0](https://doi.org/10.5285/dac20505-a56e-4beb-97ba-077eecd587c0). [BAS data documentation](https://data.bas.ac.uk/full-record.php?id=GB/NERC/BAS/PDC/01329). This example contains 17,566 epochs from 2005-03-01 to 2005-05-01 (exclusive). The project's existing along-flow projection uses an azimuth of 152.516413° clockwise from true north. The original GPS documentation is provided in `data/original_data_readme.txt`.

Times follow the project's `datetime(year,1,1) + days(doy-1)` convention. The source documentation does not specify whether the time standard is UTC or GPS; the revised version provisionally treats it as UTC. The default reproduction of the original figure instead uses the empirical time mapping described above.

- The Rutford observations are used under the [Open Government Licence v3.0](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/) listed by BAS. The supplied inputs are processed derivatives.
- The CATS2008 constants are used under the [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) listed by USAP-DC.
- The three functions in `third_party/tmd/` required by the revised version come from [TMD3](https://github.com/chadagreene/Tide-Model-Driver). The original code and MIT licence file `LICENSE` are retained; cite [Greene et al. (2024)](https://doi.org/10.21105/joss.06018). The default `RUN.m` does not depend on these helper functions.
