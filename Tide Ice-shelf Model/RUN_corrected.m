% RUN_CORRECTED Alternative physical model, with extended simulation coverage.
% Requires only MATLAB (tested with R2024b). Inputs and outputs use this folder.
% Edit physical parameters in run_example below; no parameter fitting is run.
[result, comparison, settings] = run_example(fileparts(mfilename('fullpath')));

function [result, comparison, settings] = run_example(root)
old_path = path;
restore_path = onCleanup(@() path(old_path));
addpath(fullfile(root,'core'),fullfile(root,'third_party','tmd'),'-begin');

% Rutford physical settings from the project's existing four-tide example.
settings.parameters = [1.54, 75e3, 55e3, 160e3]; % alpha, stresses (Pa)
settings.length_scale_m = 20e3;
settings.E_Pa = 9e9;
settings.A_Glen = 7e-25; % Pa^-3 s^-1; Glen exponent n=3
settings.tidal_amp = 1; % dimensionless normalized forcing multiplier
settings.dt_seconds = 300;
settings.simulation_shift_days = -12; % display only: negative = left, 0 = unshifted
settings.reference_lat_lon = [-78.8, -82.0];
settings.gps_lat_lon = [-78.141646342, -88.91234181];

data = readtable(fullfile(root,'data','Rutford_observations.csv'), ...
    'TextType','string');
time = datetime(data.time,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
tides = readtable(fullfile(root,'data','Rutford_CATS2008.csv'), ...
    'TextType','string');
assert(isequal(tides.constituent,["M2";"S2";"K1";"O1"]));
assert(all(diff(time)>seconds(0)) && height(data)>2);
assert(numel(unique(data.segment_id))==1,'Example must not cross relocations.');

% Extend the computation, rather than truncating the curve after a left shift.
model_start = min(time(1),datetime(2005,3,1)-days(settings.simulation_shift_days));
model_end = max(time(end),datetime(2005,5,1)-days(settings.simulation_shift_days));
model_time = unique([time; (model_start:seconds(settings.dt_seconds): ...
    model_end+seconds(settings.dt_seconds))']);
options = struct('Time',model_time,'PhaseConvention','greenwich_lag', ...
    'Constituents',tides.constituent,'LengthScale_m',settings.length_scale_m, ...
    'E_Pa',settings.E_Pa,'A_Glen',settings.A_Glen);
% Values of observations are deliberately not passed to the forward model.
result = Drygalski_forward_model_corrected(zeros(size(model_time)), ...
    settings.dt_seconds,tides.period_hours*3600,tides.amplitude_m, ...
    settings.tidal_amp,tides.phase_lag_deg,settings.parameters,options);
unshifted_m = interp1(seconds(model_time-model_start),result.displacement_detrend, ...
    seconds(time-model_start));
comparison = table(time,data.observed_m,unshifted_m, ...
    'VariableNames',{'time','observed_m','simulated_m'});
comparison.simulated_plot_time = time + days(settings.simulation_shift_days);
simulation = table(model_time,model_time+days(settings.simulation_shift_days), ...
    result.displacement_detrend,'VariableNames',{'model_time','display_time','simulated_m'});

fig = figure('Color','w','Position',[100 100 1150 480]);
ax = axes(fig);
% Insert plot-only NaNs so lines never bridge a missing observation interval.
gap = [false; seconds(diff(time))>1.5*settings.dt_seconds];
plot_time = NaT(height(data)+nnz(gap),1);
plot_y = nan(numel(plot_time),2);
positions = (1:height(data))' + cumsum(gap);
plot_time(positions) = time;
plot_y(positions,:) = 100*[comparison.observed_m,comparison.simulated_m];
plot(ax,plot_time,plot_y(:,1),'Color',[0 114 178]/255,'LineWidth',0.7);
hold(ax,'on');
plot(ax,simulation.display_time,100*simulation.simulated_m, ...
    'Color',[213 94 0]/255,'LineWidth',1);
model_label = sprintf('Model simulation (time shift: %+g days)',settings.simulation_shift_days);
legend(ax,{'GNSS observation',model_label},'Location','best','Box','off');
xlabel(ax,'Date (2005)'); ylabel(ax,'Along-flow displacement (cm)');
title(ax,'Rutford | M2 + S2 + K1 + O1 forcing');
xlim(ax,[datetime(2005,3,1),datetime(2005,5,1)]);
xtickformat(ax,'dd/MM'); ax.XAxis.SecondaryLabel.String = '';
xlabel(ax,'Date (dd/mm, 2005)'); grid(ax,'on'); box(ax,'on');
set(ax,'FontSize',12,'TickDir','out');

out = fullfile(root,'output','corrected');
if ~isfolder(out), mkdir(out); end
writetable(comparison,fullfile(out,'Rutford_comparison.csv'));
writetable(simulation,fullfile(out,'Rutford_simulation.csv'));
save(fullfile(out,'Rutford_result.mat'),'result','comparison','settings','tides','simulation');
exportgraphics(fig,fullfile(out,'Rutford_comparison.png'),'Resolution',200);
fprintf('Rutford: %d epochs, %s to %s.\n',height(data),string(time(1)),string(time(end)));
fprintf('Saved figure, comparison CSV and model MAT to %s\n',out);
end
