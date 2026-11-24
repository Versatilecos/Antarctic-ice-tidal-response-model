function result=Drygalski_forward_model_corrected(raw_up,dt,w_tides,A_tides,tidal_amp,phase_degrees,parameters,options)
% Deterministic physical forward model. NEVER uses raw_up values, only length.
% Original seven positional inputs retained; optional eighth options struct.
% parameters=[alpha,sigma_b0_Pa,sigma_h0_Pa,sigma_other_Pa,lag_hours].
% Fifth parameter optional, default zero. w_tides: periods in seconds.
% tidal_amp: dimensionless multiplier [0,1] of normalized ocean tide.
% options.Time: N increasing datetime output epochs; dt is integration step.
% options.PhaseConvention: 'sine_initial' (default) or 'greenwich_lag'.
% Greenwich mode requires Time and options.Constituents and uses TMD3 factors.
% options.LengthScale_m, E_Pa, A_Glen optionally override physical constants.
% No empirical harmonics, fitting, 300-day extension or arbitrary time shift.
if nargin<5||isempty(tidal_amp),tidal_amp=1;end
if nargin<6||isempty(phase_degrees),phase_degrees=zeros(size(w_tides));end
if nargin<7||isempty(parameters),parameters=[2.53,5e3,30e3,20e3];end
if nargin<8,options=struct();end
P=parameters(:)';T=w_tides(:)';A=A_tides(:)';g=phase_degrees(:)';
assert(~isempty(raw_up)&&isscalar(dt)&&isfinite(dt)&&dt>0);
assert(numel(P)==4||numel(P)==5,'Use four parameters plus optional lag hours.');
assert(all(isfinite(P))&&P(1)>0&&P(2)>=0&&P(3)>0&&P(4)>=0);
lag=0;if numel(P)==5,lag=P(5);end
assert(lag>=0,'Response lag must be nonnegative.');
assert(~isempty(T)&&numel(T)==numel(A)&&numel(T)==numel(g));
assert(all(isfinite(T)&T>0)&&all(isfinite(A)&A>=0)&&sum(A)>0&&all(isfinite(g)));
assert(isscalar(tidal_amp)&&isfinite(tidal_amp)&&tidal_amp>=0&&tidal_amp<=1, ...
    'tidal_amp is dimensionless and must be between 0 and 1.');
if isfield(options,'Time')
    output_time=options.Time(:);
    assert(isdatetime(output_time)&&numel(output_time)==numel(raw_up)&&all(~isnat(output_time)));
    output_seconds=seconds(output_time-output_time(1));assert(all(diff(output_seconds)>0));epoch=output_time(1);
else
    output_seconds=(0:numel(raw_up)-1)'*dt;output_time=[];epoch=[];
end
mode=string(getoption('PhaseConvention','sine_initial'));
assert(any(mode==["sine_initial","greenwich_lag"]));
L=getoption('LengthScale_m',25e3);E=getoption('E_Pa',9e9);AGlen=getoption('A_Glen',7e-25);
assert(isscalar(L)&&isfinite(L)&&L>0&&isscalar(E)&&isfinite(E)&&E>0&&isscalar(AGlen)&&isfinite(AGlen)&&AGlen>0);
tg=(0:dt:dt*ceil(output_seconds(end)/dt))';omega=2*pi./T;
if mode=="greenwich_lag"
    assert(~isempty(epoch)&&isfield(options,'Constituents'),'Greenwich phases require Time and Constituents.');
    % TMD astronomical helpers are supplied by the standalone RUN entry point.
    names=cellstr(lower(string(options.Constituents(:))))';assert(numel(names)==numel(T));
    [~,~,phase0,reference_omega]=tmd_constit(names);
    assert(max(abs(omega-reference_omega(:)'))<1e-10,'Periods do not match constituent names.');
    pad=max(lag,getoption('NormalizationPaddingHours',72));
    tf=(-dt*ceil(pad*3600/dt):dt:tg(end))';dates=epoch+seconds(tf);
    [~,~,p,N]=tmd_astrol(datenum(dates));[u,f]=tmd_nodal(names,p,N);
    theta=((datenum(dates)-datenum(1992,1,1))*86400)*omega+phase0(:)'+u;
    tide=sum(f.*A.*cos(theta-deg2rad(g)),2);normalization=sum(max(f,[],1).*A);
    h=interp1(tf,tidal_amp*tide/normalization,tg-lag*3600,'linear');
    ocean_unlagged=interp1(tf,tide,output_seconds,'linear');
else
    phase=(tg-lag*3600)*omega+deg2rad(g);normalization=sum(A);
    h=tidal_amp*sum(A.*sin(phase),2)/normalization;
    ocean_unlagged=sum(A.*sin(output_seconds*omega+deg2rad(g)),2);
end
assert(all(isfinite(h))&&max(abs(h))<=1+1e-10);
shape=2*max(0,(h+1)/2).^P(1)-1;
sigma=P(2)*shape-P(3)*h;
eff2=.5*P(2)^2*shape.^2+.5*P(3)^2*h.^2+P(4)^2;
elastic=L/E*(sigma-sigma(1));
viscous_velocity=L*(2/3)*AGlen*eff2.*sigma;
if numel(tg)>1,viscous=cumtrapz(tg,viscous_velocity);else,viscous=0;end
displacement=sample(elastic+viscous);
if numel(tg)>1,dsigma=gradient(sigma,dt);else,dsigma=0;end
% Elastic integral above is exact; gradient is for diagnostics only.
velocity=L/E*dsigma+viscous_velocity;days_out=output_seconds/86400;
if numel(displacement)>1
    B=[ones(size(days_out)),days_out-mean(days_out)];trend=B*(B\displacement);
else,trend=displacement;end
result=struct('time',output_time,'time_seconds',output_seconds,'time_days',days_out, ...
    'h_star',sample(h),'ocean_tide_unlagged',ocean_unlagged,'normalization_amplitude',normalization, ...
    'sigma',sample(sigma),'sigma_eff',sqrt(sample(eff2)),'dsigma_dt',sample(dsigma), ...
    'velocity',sample(velocity),'displacement',displacement,'displacement_trend',trend, ...
    'displacement_detrend',displacement-trend,'displacement_detrend_cm',100*(displacement-trend), ...
    'elastic_displacement',sample(elastic),'viscous_displacement',sample(viscous), ...
    'w_tides',T,'A_tides',A,'phase_degrees',g,'tidal_amp',tidal_amp,'phase_lag_deg',NaN);
result.params=struct('alpha',P(1),'sigma_b0',P(2),'sigma_h0',P(3),'sigma_other',P(4), ...
    'beta',P(2)/P(3),'gamma',P(4)/P(3),'lag_hours',lag,'length_scale',L,'E',E,'A_Glen',AGlen,'n',3);
result.metadata=struct('version','2026-09-09 corrected deterministic source model', ...
    'phase_convention',mode,'integration_step_seconds',dt,'observed_values_used',false, ...
    'time_assumption','Input dates assumed UTC; source standard needs confirmation', ...
    'normalization_note','Amplitude ratios retained; sigma_h0 sets stress scale', ...
    'phase_lag_note','Broadband FFT phase no longer reported; use constituent diagnostics');
    function value=getoption(name,default)
        if isfield(options,name),value=options.(name);else,value=default;end
    end
    function value=sample(x)
        if numel(tg)==1,value=x(1);else,value=interp1(tg,x,output_seconds,'linear');end
    end
end
