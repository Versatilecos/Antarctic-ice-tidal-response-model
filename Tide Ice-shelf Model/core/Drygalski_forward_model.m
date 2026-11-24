function result = Drygalski_forward_model(raw_up, sample_interval_seconds, w_tides, A_tides, tidal_amp, phase_degrees,parameters)
%DRYGALSKI_FORWARD_MODEL Forward model for ice-shelf tidal response.
%
% result = Drygalski_forward_model(raw_up, sample_interval_seconds, w_tides, A_tides)
% result = Drygalski_forward_model(raw_up, sample_interval_seconds, w_tides, A_tides, tidal_amp)
% result = Drygalski_forward_model(raw_up, sample_interval_seconds, w_tides, A_tides, tidal_amp, phase_degrees)
%
% Inputs
%   raw_up                  Raw displacement series. Its values are not used
%                           by the forward model; its length defines the
%                           simulated output length.
%   sample_interval_seconds Model time step, in seconds. To exactly match
%                           Model_Drygalski3.m, use
%                           (10*60*3600)/(numel(raw_up)-1).
%   w_tides                 Tide periods, in seconds. Example:
%                           [24*3600*0.9971, 24*3600*0.5175].
%   A_tides                 Tide amplitudes or weights. Must match w_tides.
%   tidal_amp               Optional total tidal height scale, in meters.
%   phase_degrees           Optional initial phase of each constituent, in
%                           degrees. Must match w_tides. Default is zero.
%   parameters              Physical parameters, including alpha and sigma
%   choice params
%   alpha = 1.54;                   %Exponent of buttressing stress w.r.t tidal height
%   tidal_amp = 3;                  %Amplitude of combined tides
%   sigma_b0 = 75e3;                %Amplitude of tidal changes in buttressing stress
%   sigma_h0 = 55e3;                %Amplitude of tidal changes in hydrostatic stress
%   sigma_other = 160e3;            %Total background stresses
%   length_scale = 25e3;            %Length scale between two stations to calculate strain
% Output
%   result                  Struct containing simulated tide, stress,
%                           velocity, displacement, and diagnostic values.

if nargin < 5 || isempty(tidal_amp)
    tidal_amp = 1;
end

if nargin < 6 || isempty(phase_degrees)
    phase_degrees = zeros(size(w_tides));
end

raw_up = raw_up(:);
w_tides = w_tides(:).';
A_tides = A_tides(:).';
phase_degrees = phase_degrees(:).';

if isempty(raw_up)
    error('Drygalski_forward_model:EmptyInput', 'raw_up must contain at least one sample.');
end

if ~isscalar(sample_interval_seconds) || sample_interval_seconds <= 0
    error('Drygalski_forward_model:InvalidSampleInterval', ...
        'sample_interval_seconds must be a positive scalar.');
end

if isempty(w_tides) || numel(w_tides) ~= numel(A_tides)
    error('Drygalski_forward_model:InvalidTideParameters', ...
        'w_tides and A_tides must be non-empty vectors with the same length.');
end

if numel(phase_degrees) ~= numel(w_tides)
    error('Drygalski_forward_model:InvalidTidePhases', ...
        'phase_degrees must have the same length as w_tides.');
end

if any(w_tides <= 0)
    error('Drygalski_forward_model:InvalidTidePeriods', 'All tide periods must be positive.');
end

if sum(A_tides) == 0
    error('Drygalski_forward_model:InvalidTideAmplitudes', 'The sum of A_tides must not be zero.');
end

if any(phase_degrees < 0 | phase_degrees > 360)
    error('Drygalski_forward_model:InvalidTidePhases', ...
        'phase_degrees values must be between 0 and 360.');
end

phase_degrees = mod(phase_degrees, 360);
phase_radians = phase_degrees * pi / 180;

%% 控制额外反演天数 
% 总反演天数=raw_up天数+extra_days
extra_days=300;
extra_samples = round(extra_days * 24*3600 / sample_interval_seconds);

%% Physical parameters
%length_scale = 25e3;

E = 9e9;
A_Glen = 7e-25;
n = 3;

% alpha = 2.53;
% sigma_b0 = 5e3;
% sigma_h0 = 30e3;
% sigma_other = 20e3;
alpha = parameters(1);
sigma_b0 = parameters(2);
sigma_h0 = parameters(3);
sigma_other = parameters(4);
length_scale = parameters(5);

beta = sigma_b0 / sigma_h0;
gamma = sigma_other / sigma_h0;

%% Time settings
num_samples = numel(raw_up) + extra_samples;
nt = num_samples - 1;
dt = sample_interval_seconds;
tf = nt * dt;
time_seconds = (0:nt).' * dt;
time_days = time_seconds / (24 * 3600);

%% Integrate to a quasi-steady initial state
disp_value = 0;
spinup_dt = 3600 * 24 / 4;
spinup_time = spinup_dt * 5e3;

for t = 0:spinup_dt:spinup_time
    [~, sigma, dsigma_dt, sigma_eff] = local_tide_stress( ...
        t, w_tides, A_tides, phase_radians, tidal_amp, alpha, beta, gamma, sigma_h0);

    velocity = length_scale * ((1 / E) .* dsigma_dt + ...
        (2 / 3) * A_Glen .* (sigma_eff.^(n - 1)) .* sigma);
    disp_value = disp_value + velocity * dt;
end

%% Integrate over the requested data length
hs = zeros(num_samples, 1);
sigmas = zeros(num_samples, 1);
disps = zeros(num_samples, 1);
vels = zeros(num_samples, 1);
sigma_effs = zeros(num_samples, 1);
dsigma_dts = zeros(num_samples, 1);

for k = 1:num_samples
    t = time_seconds(k);
    [h_star, sigma, dsigma_dt, sigma_eff] = local_tide_stress( ...
        t, w_tides, A_tides, phase_radians, tidal_amp, alpha, beta, gamma, sigma_h0);

    velocity = length_scale * ((1 / E) .* dsigma_dt + ...
        (2 / 3) * A_Glen .* (sigma_eff.^(n - 1)) .* sigma);

    disp_value = disp_value + velocity * dt;

    hs(k) = h_star;
    sigmas(k) = sigma;
    disps(k) = disp_value;
    vels(k) = velocity;
    sigma_effs(k) = sigma_eff;
    dsigma_dts(k) = dsigma_dt;
end

if num_samples > 1
    trend = polyval(polyfit(time_days, disps, 1), time_days);
else
    trend = disps;
end
disps_detrend = disps - trend;

result = struct();
result.time_seconds = time_seconds;
result.time_days = time_days;
result.h_star = hs;
result.sigma = sigmas;
result.sigma_eff = sigma_effs;
result.dsigma_dt = dsigma_dts;
result.velocity = vels;
result.displacement = disps;
result.displacement_detrend = disps_detrend;
result.displacement_detrend_cm = 100 * disps_detrend;
result.w_tides = w_tides;
result.A_tides = A_tides;
result.phase_degrees = phase_degrees;
result.tidal_amp = tidal_amp;
result.params = struct( ...
    'length_scale', length_scale, ...
    'E', E, ...
    'A_Glen', A_Glen, ...
    'n', n, ...
    'alpha', alpha, ...
    'sigma_b0', sigma_b0, ...
    'sigma_h0', sigma_h0, ...
    'sigma_other', sigma_other, ...
    'beta', beta, ...
    'gamma', gamma);

result.phase_lag_deg = local_phase_lag(hs, disps_detrend);

end

function [h_star, sigma, dsigma_dt, sigma_eff] = local_tide_stress( ...
    t, w_tides, A_tides, phase_radians, tidal_amp, alpha, beta, gamma, sigma_h0)

phase = 2 * pi * t ./ w_tides + phase_radians;
amp_sum = sum(A_tides);

h_star_unit = sum(A_tides .* sin(phase)) / amp_sum;
h_star = h_star_unit * tidal_amp;

% This follows Model_Drygalski3.m: the stress derivative uses an unweighted
% average of component frequency terms, not the exact derivative of h_star.
dhs_dt = sum((2 * pi ./ w_tides) .* cos(phase)) / numel(w_tides);

normalized_tide = 0.5 * (h_star + 1);
buttressing_shape = 2 * normalized_tide.^alpha - 1;

sigma = sigma_h0 * (beta * buttressing_shape - h_star);
dsigma_dt = sigma_h0 .* dhs_dt .* ...
    (alpha * beta * normalized_tide.^(alpha - 1) - 1);

sigma_eff = sigma_h0 .* sqrt( ...
    0.5 * beta.^2 * buttressing_shape.^2 + ...
    0.5 * h_star.^2 + ...
    gamma.^2);

end

function phase_lag_deg = local_phase_lag(h_star, displacement_detrend)

if numel(h_star) < 2
    phase_lag_deg = NaN;
    return
end

L = length(h_star);
NFFT = 2^nextpow2(L);
fft_h = fft(h_star, NFFT) / L;
fft_d = fft(displacement_detrend, NFFT) / L;
PSD = 2 * abs(fft_h(1:NFFT / 2 + 1));
[~, idx_primary_tide] = max(PSD);

if abs(fft_d(idx_primary_tide)) == 0
    phase_lag_deg = NaN;
else
    phase_lag_deg = 180 * angle(fft_h(idx_primary_tide) / fft_d(idx_primary_tide)) / pi;
end

end
