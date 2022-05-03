
clearvars;
%%%%%%%%%%%%%%%%%% Simulation 5 - Cylindrical Cortical Bone Slab %%%%%%%%%%%%
% =========================================================================
% DEFINE LITERALS
% =========================================================================
    

% medium parameters
c0              = 1500;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]
u               = 0.04;
cbrain          = 1560;     % speed of sound [m/s]
rhobrain        = 1040; 
cskullcortical  = 2800;     % speed of sound [m/s]
rhoskullcortical= 1850; 
%% attenuation coefficients
alphabrain = 1.2;
alphatrab = 32;
alphacort = 16;
alpha0 = 0;


% source parameters
source_f0       = 500000;   % source frequency [Hz]
source_roc      = 15e-3;    % bowl radius of curvature [m]
source_diameter = 10e-3;    % bowl aperture diameter [m]
source_mag      = rho0*c0*u;% source pressure [Pa]
offset          = source_roc - sqrt((source_roc^2)-(source_diameter/2)^2);

% grid parameters
axial_size      = 64e-3;    % total grid size in the axial dimension [m]
lateral_size    = 40e-3;   % total grid size in the lateral dimension [m]

ppw             = 12;        % number of points per wavelength
t_end           = 120e-6;    % total compute time [s] (this must be long enough to reach steady state)
record_periods  = 3;        % number of periods to record
cfl             = 0.05;     % CFL number

% =========================================================================
% RUN SIMULATION
% =========================================================================

% --------------------
% GRID
% --------------------

% calculate the grid spacing based on the PPW and F0
dx = c0 / (ppw * source_f0);   % [m]
%dx =10e-04
%dx = 0.00027
source_x_offset = round(offset/dx); %grid points for offset of source 

% compute the size of the grid
Nx = roundEven(axial_size / dx) 
%+ source_x_offset;
Ny = roundEven(lateral_size / dx);

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx);

% compute points per temporal period
PPP = round(ppw / cfl);

% compute corresponding time spacing
dt = 1 / (PPP * source_f0);

% create the time array using an integer number of points per period
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% calculate the actual CFL and PPW
disp(['PPW = ' num2str(c0 / (dx * source_f0))]);
disp(['CFL = ' num2str(c0 * dt / dx)]);

% --------------------
% SOURCE   - Planar Piston
% --------------------

% create time varying source
source_sig = createCWSignals(kgrid.t_array, source_f0, source_mag, 0);

%create empty kWaveArray
karray = kWaveArray('Axisymmetric', true, 'BLITolerance', 0.01);
source_diam     = 0.02; 
% add line shaped element
karray.addLineElement([kgrid.x_vec(1), -source_diam/2], [kgrid.x_vec(1), source_diam/2]);


% assign binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% assign source signals
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

% plot grid weights
grid_weights = karray.getArrayGridWeights(kgrid);
figure;
imagesc(kgrid.y_vec - kgrid.y_vec(1), kgrid.x_vec - kgrid.x_vec(1), grid_weights);
axis image;
title('Off-grid source weights');
%%


% --------------------
% MEDIUM
% --------------------

%Create geometry
medskull = zeros(Nx,Ny);
medskull(120:146,1:140) = 1;
imagesc(medskull);
axis equal
% 
% % 
% Define Medium parameters

medium.sound_speed= cskullcortical*medskull;
medium.sound_speed(medium.sound_speed == 0) = c0;

medium.density = rhoskullcortical*medskull;
medium.density(medium.density == 0) = rho0;

medium.alpha_power = 2;

medium.alpha_coeff = alphacort*medskull;
medium.alpha_coeff(medium.alpha_coeff == 0) = alpha0;

%% Dimensions of geometry in mm %%
cuboidstart = dx*120*1e3
cuboidend = dx* 146*1e3
cuboidwidth = dx*26*1e3
cuboiddiam = dx*140*1e3
%%
% --------------------
% SENSOR
% --------------------

% set sensor mask to record central plane, not including the source point
sensor.mask = zeros(Nx, Ny);
sensor.mask(:, :) = 1;

% record the pressure
sensor.record = {'p'};

% --------------------
% SIMULATION
% --------------------

% set input options
input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off'};

% run code

sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor, ...
    input_args{:}, ...
    'DataCast', 'single', ...
    'PlotScale', [-1, 1] * source_mag);
%%

% extract amplitude from the sensor data
amp = extractAmpPhase(sensor_data.p, 1/kgrid.dt, source_f0, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);
%%

% reshape data
amp = reshape(amp, Nx, Ny);
%- source_x_offset - 1
% extract pressure on axis
amp_on_axis = amp(:, 1);

% define axis vectors for plotting
x_vecc = kgrid.x_vec(source_x_offset + 2:end, :) - kgrid.x_vec(source_x_offset + 1);
y_vec = kgrid.y_vec - kgrid.y_vec(1);
%%

% =========================================================================
% Extracting Pressure and Time Data
% =========================================================================
x_vec = kgrid.x_vec(2:end, :) - kgrid.x_vec(1);

timeperiod= 0:dt:t_end-dt;
periodofwave = 1/source_f0;
cyclelength = periodofwave/dt;
timeindex = Nt-(cyclelength-1);
wavetime = timeperiod(timeindex:Nt);
samplingfrequency = 1/dt;
freqv = (0:cyclelength/2-1)*samplingfrequency/cyclelength;
FFTarray = [];
FFTmag = [];
for i=1:Nx-1
    pressuretime = sensor_data.p(i, 1:Nt);
    cycleextract = pressuretime(timeindex:Nt);
    FFT = fft(cycleextract);

    %symmetric so remove second half
    FFT = FFT(1:(round(cyclelength,0))/2);
    %magnitude of fft
    mFFT = (abs(FFT));
    scaledFFT = mFFT*(2/PPP);
    [FFTpeak, FFTindex] = max(scaledFFT);
    peakfreq = freqv(FFTindex);
    FFTarray = [FFTarray,peakfreq];
    FFTmag = [FFTmag,FFTpeak];
end
% =========================================================================
% VISUALISATION
% =========================================================================

% plot fft 
figure;
plot(x_vec*10^3,FFTmag/10^3)
ylabel('Magnitude of Acoustic Pressure at Peak Frequency P(f) [kPa]')
xlabel('Axial Position (mm)')
title('Acoustic Pressure at Peak Frequency on Axial Positions after reaching steady-state with a bowl piston')
%%


% % =========================================================================
% % VISUALISATION
% % =========================================================================

% duplicate the pressure field
amp = [flip(amp(:, 2:end), 2), amp];
y_vec = [-flip(y_vec(2:end)); y_vec];

% plot the pressure field 
figure;
imagesc(1e3 * x_vecc, 1e3 * y_vec, 1e-3*amp.');
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
axis image;
title('Acoustic Pressure Field with a spherical heterogeneity and focused bowl source generated using k-Wave');
a = colorbar;
a.Label.String = 'Acoustic Pressure (kPa)';

%%
%%%% Import solution from alternative solver OptimUS
load('optimus5.mat','p_amp','p_phase')
y_vec = kgrid.y_vec - kgrid.y_vec(1);
y_vec = [-flip(y_vec(2:end)); y_vec];

optx = p_amp;
surf(optx)
figure;

dxx = 5e-4;
optxx = optx(:,71);

xarrayy = linspace(0,((234-1)*dxx),193);

xarray = linspace(-dxx/2,((241-0.5)*dxx),241);
yarray = linspace(((-96)*dxx),((96)*dxx),141);
imagesc(xarray*10^3,yarray*10^3,1e-3*optx');
ylabel('Lateral Position [mm]')
xlabel('Axial Position [mm]')
title('Optimus pressure field')
colorbar
axis equal
figure;
plot(xarray*10^3,optxx/10^3);
hold on 
plot((x_vec*10^3),FFTmag/10^3)
ylabel('Magnitude of Acoustic Pressure at Peak Frequency P(f) [kPa]')
xlabel('Axial Position (mm)')
title('On-Axis Acoustic Pressure with a planar piston and cylindrical skull heterogeneity')

hold on
%%%% Import Alternative Solver Solution - k-Wave Benchmark
load('benchmark.mat','p_amp','p_phase')
y_vec = kgrid.y_vec - kgrid.y_vec(1);
y_vec = [-flip(y_vec(2:end)); y_vec];

optx1 = p_amp;
dxx = 5e-4;
benchh = optx1(:,71);
plot(xarray*10^3,benchh/10^3);
xlim([0 64])
legend('Optimus Solution','Calculated k-Wave Solution', 'Benchmark k-Wave Solution');
hold off
surf(optx)
optxxx = optx(1:129,12:130)
xarray(129)
yarray(130)
yarray(12)
imagesc(xarray(1:129)*10^3,yarray(12:130)*10^3,1e-3*optxxx');
xlim([0.3 64])
axis equal;
a = colorbar;
a.Label.String = 'Acoustic Pressure (kPa)';
%%
%%% Error Analysis
figure
vq2 = interp1(xarray,optxx,x_vec(2:end));
plot((x_vec(2:end)*10^3),vq2/10^3)
hold on 
plot((x_vec(2:end)*10^3),FFTmag(2:end)/10^3)
ylabel('Magnitude of Acoustic Pressure at Peak Frequency P(f) [kPa]')
xlabel('Axial Position (mm)')
title('Acoustic Pressure at Peak Frequency on Axial Positions after reaching steady-state with a bowl piston')
legend('optimus','k-wave')
xlim([0 59.5])
hold off
pctError = mean((abs(transpose(FFTmag(2:end))-vq2)./vq2)*100)
D = (vq2-(transpose(FFTmag(2:end)))).^2;
S1= sum(D);
S2 = sum((vq2).^2);
ctest = 100*(sqrt(S1/S2))
standarddev = std((abs(transpose(FFTmag(20:250))-vq2(19:249))./vq2(19:249))*100)