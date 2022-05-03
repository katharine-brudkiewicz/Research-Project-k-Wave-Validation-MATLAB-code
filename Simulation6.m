
clearvars;

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
axial_size      = 60e-3;    % total grid size in the axial dimension [m]
lateral_size    = 25e-3;
% computational parameters
ppw             = 11;        % number of points per wavelength
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
% SOURCE Planar Piston
% --------------------

% create time varying source
source_sig = createCWSignals(kgrid.t_array, source_f0, source_mag, 0);
% 

%create empty kWaveArray
karray = kWaveArray('Axisymmetric', true, 'BLITolerance', 0.01);
source_diam     = 0.02; 
%add line shaped element
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

% 
% --------------------
% MEDIUM
% --------------------
%%% Code to create concentric sphere geometry %%%%%%
nsphere = Ny;
spherestart = 36;
Psphere = [randi(1,nsphere,1); zeros(nsphere, 1)];
[Xsphere, Ysphere] = meshgrid(1:nsphere, flip(1:nsphere));
rsphere = round(sqrt(Xsphere.^2 + Ysphere.^2));
Asphere = reshape(Psphere(rsphere(:)), nsphere, nsphere);
Bsphere = flip(Asphere);
fillmatrix = (Nx-(length(Asphere)*2))-spherestart;
SC = [zeros(spherestart,Ny);Asphere;Bsphere;zeros(fillmatrix,Ny)];

nsphere1 = Ny-2;
spherestart1 = 38;
Psphere1 = [randi(1,nsphere1,1); zeros(nsphere1, 1)];
[Xsphere, Ysphere] = meshgrid(1:nsphere1, flip(1:nsphere1));
rsphere = round(sqrt(Xsphere.^2 + Ysphere.^2));
Csphere = reshape(Psphere1(rsphere(:)), nsphere1, nsphere1);
Dsphere = flip(Csphere);
fillmatrix1 = (Nx-(length(Csphere)*2))-spherestart1;
SC1 = [zeros(spherestart1,nsphere1);Csphere;Dsphere;zeros(fillmatrix1,nsphere1)];
SC1 = [SC1,zeros(Nx,2)];

SCO = SC-SC1;
SCO(SCO==1)=2;
SCO(SCO==0)=1;
SCO(SCO==2)=0;


% define the medium properties
medium.sound_speed = cbrain*SC; % [m/s]
medium.sound_speed(medium.sound_speed == 0) = c0; % [m/s]
medium.sound_speed = medium.sound_speed.*SCO;
medium.sound_speed(medium.sound_speed == 0) = cskullcortical;
imagesc(medium.sound_speed)

medium.density = rhobrain*SC; % [kg/m^3]
medium.density(medium.density == 0) = rho0; % [kg/m^3]
medium.density = medium.density.*SCO; 
medium.density(medium.density == 0) = rhoskullcortical;

diameterofsphere = length(Asphere)*2*dx*10^3 % [m]
locationofspherestart = spherestart*dx*10^3
locationofsphereend = (diameterofsphere)+locationofspherestart
centreofsphere = locationofspherestart + (0.5*diameterofsphere)
diameterofsphere1 = length(Csphere)*2*dx*10^3 % [m]
locationofspherestart1 = spherestart1*dx*10^3
locationofsphereend1 = (diameterofsphere1)+locationofspherestart1
centreofsphere1 = locationofspherestart1 + (0.5*diameterofsphere1)

% 
medium.alpha_power = 2;
medium.alpha_coeff = medium.sound_speed;
medium.alpha_coeff(medium.alpha_coeff == c0) = alpha0;
medium.alpha_coeff(medium.alpha_coeff == cskullcortical) = alphacort;
medium.alpha_coeff(medium.alpha_coeff == cbrain) = alphabrain;

imagesc(medium.alpha_coeff)
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
% Extracting Pressure and Time Data using FFT analysis
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
caxis([0 160])
%%
%%%%%%%%% Importing Alternative Solver Solution OptimUS %%%%%
load('optimus2.mat','p_amp','p_phase')
y_vec = kgrid.y_vec - kgrid.y_vec(1);
y_vec = [-flip(y_vec(2:end)); y_vec];
optx = reshape(p_amp,234,183);

surf(optx)
figure;
dxx = 2.727272727272727e-04;

optxx = optx(:,92);

xarrayy = linspace(0,((234-1)*dxx),193);

xarray = linspace(-dxx/2,((234-0.5)*dxx),234);
yarray = linspace(((-96)*dxx),((96)*dxx),141);
imagesc(xarray*10^3,yarray*10^3,1e-3*optx');
ylabel('Lateral Position [mm]')
xlabel('Axial Position [mm]')
title('Optimus pressure field')
a = colorbar;
a.Label.String = 'Acoustic Pressure (kPa)';

axis equal
%%
xarray(end)
figure;
plot(xarray*10^3,optxx/10^3);
hold on 
plot((x_vec*10^3),FFTmag/10^3)
ylabel('Magnitude of Acoustic Pressure at Peak Frequency P(f) [kPa]')
xlabel('Axial Position (mm)')
title('On-Axis Acoustic Pressure with a planar piston and cylindrical skull heterogeneity')


xlim([0 60])
legend('Optimus Solution','Calculated k-Wave Solution');
hold off
surf(optx)
%%
%% Visualisation of Optmius
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
%%% Calculating Error of Solution
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
standarddev = std((abs(transpose(FFTmag(2:end))-vq2)./vq2)*100)
figure
pctError = (abs(transpose(FFTmag)-optxx(2:end))./optxx(2:end))*100;
%pctError = mean(100 * abs(FFTmag-vq2) ./ vq2);
plot((x_vec*10^3),pctError)

ylabel('Percentage error of acoustic pressure')
xlabel('Axial Position (mm)')
title('Acoustic Pressure at Peak Frequency on Axial Positions after reaching steady-state with a bowl piston')

xlim([0 59.5])
ppppp = pctError(1:50)
percentageerror = mean(ppppp)