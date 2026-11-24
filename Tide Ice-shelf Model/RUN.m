% RUN Reproduce the original Rutford M2+S2+K1+O1 figure. Click Run.
[result, comparison, settings, simulation] = run_example(fileparts(mfilename('fullpath')));

function [result, comparison, settings, simulation] = run_example(root)
old_path = path; restore_path = onCleanup(@() path(old_path));
addpath(fullfile(root,'core'),'-begin');
settings.parameters = [1.54,75e3,55e3,160e3,20e3]; % alpha, three stresses (Pa), L (m)
settings.dt_seconds = 300;
settings.tidal_amp = 1;
settings.simulation_shift_days = -4.3; % original figure; negative = left
settings.start_time = datetime(2005,3,1);
settings.end_time = datetime(2005,5,1);
settings.extra_days = 300; % retained inside the original core
settings.reference_lat_lon = [-78.8,-82.0];
settings.gps_lat_lon = [-78.141646342,-88.91234181];
settings.model_version = 'Original Rutford figure reproduction';
settings.phase_convention = 'Original script: CATS lag values used as sine initial phases';

S = load(fullfile(root,'data','Rutford_original_observations.mat'));
time = S.time(:); obs_cm = S.observed_cm(:);
tides = readtable(fullfile(root,'data','Rutford_CATS2008.csv'),'TextType','string');
assert(isequal(tides.constituent,["M2";"S2";"K1";"O1"]));
% Match the original script exactly, including amplitude units and phases.
% Observation values are not used by the model; only the record length is used.
result = Drygalski_forward_model(zeros(size(time)),settings.dt_seconds, ...
    3600*tides.period_hours,100*tides.amplitude_m,settings.tidal_amp, ...
    tides.phase_lag_deg,settings.parameters);
model_time = settings.start_time + days(result.time_days);
plot_time = model_time + days(settings.simulation_shift_days);
selected = plot_time>=settings.start_time & plot_time<=settings.end_time;
assert(any(selected));
assert(plot_time(1)<=settings.start_time && plot_time(end)>=settings.end_time, ...
    'The simulation does not cover the display window after shifting.');
settings.baseline_cm = mean(100*result.displacement_detrend(selected));
sim_cm = 100*result.displacement_detrend-settings.baseline_cm;
simulation = table(model_time,plot_time,sim_cm/100, ...
    'VariableNames',{'model_time','display_time','simulated_m'});
aligned_m = interp1(seconds(plot_time-settings.start_time),sim_cm/100, ...
    seconds(time-settings.start_time),'linear');
comparison = table(time,obs_cm/100,aligned_m, ...
    'VariableNames',{'time','observed_m','simulated_shifted_m'});

fig = figure('Color','w','Units','centimeters','Position',[2 2 18.2 7.5]);
ax = axes(fig,'Position',[.12,.25,.85,.60]);
plot(ax,time,obs_cm,'Color',[35 83 112]/255,'LineWidth',.55); hold(ax,'on');
plot(ax,plot_time(selected),sim_cm(selected),'Color',[191 126 65]/255,'LineWidth',1);
xlim(ax,[settings.start_time,settings.end_time]); ylim(ax,[-30 30]);
yticks(ax,[-30 -15 0 15 30]);
xticks(ax,[datetime(2005,3,[1 15 29]),datetime(2005,4,[12 26])]);
xticklabels(ax,{'01 Mar','15 Mar','29 Mar','12 Apr','26 Apr'});
ax.XAxis.SecondaryLabel.Visible = 'off';
ylabel(ax,'Displacement (cm)'); xlabel(ax,'2005');
text(ax,0,1.10,'Rutford   |   M_2 + S_2 + K_1 + O_1', ...
    'Units','normalized','FontSize',9,'FontWeight','bold');
set(ax,'FontName','Arial','FontSize',8,'TickDir','out','Box','off','LineWidth',.6);
lg = legend(ax,{'GPS observation','Model simulation'}, ...
    'Orientation','horizontal','Box','off');
lg.Units = 'normalized'; lg.Position = [.28 .025 .55 .05];
ax.Position = [.12,.25,.85,.60];
fig.UserData = settings;

out = fullfile(root,'output'); if ~isfolder(out), mkdir(out); end
writetable(comparison,fullfile(out,'Rutford_comparison.csv'));
writetable(simulation,fullfile(out,'Rutford_simulation.csv'));
save(fullfile(out,'Rutford_result.mat'),'result','comparison','settings','simulation','tides');
exportgraphics(fig,fullfile(out,'Rutford_comparison.png'),'Resolution',600);
exportgraphics(fig,fullfile(out,'Rutford_comparison.pdf'),'ContentType','vector');
savefig(fig,fullfile(out,'Rutford_comparison.fig'));
fprintf('Original Rutford reproduction: %.2f simulation days, shift %+g days.\n', ...
    result.time_days(end),settings.simulation_shift_days);
fprintf('Full display window: %s to %s. Saved to %s\n', ...
    string(settings.start_time),string(settings.end_time),out);
end
