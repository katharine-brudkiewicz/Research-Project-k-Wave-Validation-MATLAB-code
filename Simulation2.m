
%clearvars;
%%%%%%%%%%%%%%%%  SIMULATION 2 - FOCUSED BOWL HOMOGENEOUS NO ATTENUATION %%
% =========================================================================
% DEFINE LITERALS
% =========================================================================
    

% medium parameters
c0              = 1500;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]
u               = 0.04;

% source parameters
source_f0       = 500000;   % source frequency [Hz]
source_roc      = 15e-3;    % bowl radius of curvature [m]
source_diameter = 20e-3;    % bowl aperture diameter [m]
source_mag      = rho0*c0*u;% source pressure [Pa]
offset          = source_roc - sqrt((source_roc^2)-(source_diameter/2)^2);

% grid parameters
axial_size      = 64e-3;    % total grid size in the axial dimension [m]
lateral_size    = 32e-3;  % total grid size in the lateral dimension [m]

% computational parameters
ppw             = 13;        % number of points per wavelength
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



% compute the size of the grid
Nx = roundEven(axial_size / dx) 

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
% SOURCE
% --------------------

% create time varying source
source_sig = createCWSignals(kgrid.t_array, source_f0, source_mag, 0);

% create empty kWaveArray
karray = kWaveArray('Axisymmetric', true, 'BLITolerance', 0.01, 'UpsamplingRate', 10);

% set bowl position and orientation
arc_pos = [kgrid.x_vec(1), 0];
focus_pos = [kgrid.x_vec(end), 0];


% add arc shaped element
karray.addArcElement(arc_pos, source_roc, source_diameter, focus_pos);

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
    
% --------------------
% MEDIUM
% --------------------

% assign medium properties
medium.sound_speed = c0;
medium.density = rho0;

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

% extract amplitude from the sensor data
amp = extractAmpPhase(sensor_data.p, 1/kgrid.dt, source_f0, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% reshape data
amp = reshape(amp, Nx, Ny);

% extract pressure on axis
amp_on_axis = amp(:, 1);

% define axis vectors for plotting
x_vec = kgrid.x_vec(source_x_offset + 2:end, :) - kgrid.x_vec(source_x_offset + 1);
y_vec = kgrid.y_vec - kgrid.y_vec(1);


% % =========================================================================
% % Extracting Pressure and Time Data  - FFT ANALYSIS
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
    FFT = FFT(1:cyclelength/2);
    
    %magnitude of fft
    mFFT = (abs(FFT));
    scaledFFT = mFFT*(2/PPP);
    [FFTpeak, FFTindex] = max(scaledFFT);
    peakfreq = freqv(FFTindex);
    FFTarray = [FFTarray,peakfreq];
    FFTmag = [FFTmag,FFTpeak];
end



%%%%%%-----------%%
%% Visualisation%%%%

% plot the pressure field 

% duplicate the pressure field
amp = [flip(amp(:, 2:end), 2), amp];
y_vec = [-flip(y_vec(2:end)); y_vec];
figure;
imagesc(1e3 * x_vec, 1e3 * y_vec, 1e-3*amp.');
%xlabel('Axial Position [mm]');
%ylabel('Lateral Position [mm]');
axis equal;
%title('Acoustic Pressure Field with a heterogeneous medium of complex cranial geometry and focused bowl source generated using k-Wave');
a = colorbar;
a.Label.String = 'Acoustic Pressure (kPa)';
%%
%%%==================================%%%
%$% Rayleigh's Integral Discretised %%$$$%


%% Variables %%

w = 311; % defines computational grid size
c = 1500;%speed of sound
f = 500000; %frequency

omega = (2*pi*f);
k = omega/c;% wave number
ww = c/f; %wavelength = speed/frequency
p0 = 1000; %density of water

%% 3D coordinates %%
u = [0.04,0,0];%velocity vector
n = [1, 0, 0]; %normal vector with positive x axis being the axis of propagation
v = dot(u,n); %dot product of velocity vector with normal vector

x2=linspace(0,0.064,w);
y2=linspace(-0.032,0.032,w);
z2=linspace(0,0.064,w);
[X,Y,Z]=meshgrid(x2,y2,z2);
sz = size(X);
P = zeros(sz); %creates matrix for receiver coordinates to sum all the values calculated below for each coordinate
%% Focused Bowl Source %%
[XA,YA,ZA] = sphere(300); %resolution of transducer defined here
radiusofcurvature = 0.015;
diameterofaperture = 0.02;
polarangle = asin((diameterofaperture/2)/radiusofcurvature);
h = radiusofcurvature - sqrt((radiusofcurvature^2)-(diameterofaperture/2)^2);
bowlarea = pi*((diameterofaperture/2)^2+h^2);
sf=surf(radiusofcurvature.*XA,radiusofcurvature.*YA,radiusofcurvature.*ZA);
Xa = get(sf, 'XData')+radiusofcurvature;
Ya = get(sf, 'YData');
Za = get(sf, 'ZData');
axis equal
spherel = length(Xa)^2;
xa=reshape(Xa,spherel,1);
ya=reshape(Ya,spherel,1);
za=reshape(Za,spherel,1);
q = length(xa);
spherematrix = [xa,ya,za];
indices = find(spherematrix(:,1)>h);
spherematrix(indices,:)=[];
xb = spherematrix(:,1);
yb = spherematrix(:,2);
zb = spherematrix(:,3);
plot3(xb,yb,zb,'o','MarkerSize',1)
axis equal

%% Pressure calculation without attenuation %%
g = bowlarea/length(xb);
t = 1;

for j=1:length(yb)
        P = P + (((1i*p0*c*k)/(2*pi)).*((((exp((sqrt((((xb(j)-X).^2)+(((yb(j))-Y).^2)+((zb(j))-Z).^2))).*k.*-1i))./(sqrt(((xb(j)-X).^2)+(((yb(j))-Y).^2)+((zb(j)-Z).^2)))).*v.*g)));
end

PP = abs(P);
%%%% Visualisation
pd = squeeze(PP(156,:,1));


plot(x2*10^3,pd/10^3);
title('Comparison of the Magnitude of Acoustic Pressure of Bowl Source from Rayleigh Integral and k-Wave')
xlabel('Axial Position (mm)')
ylabel('Acoustic Pressure [kPa]')
hold on
plot(x_vec(1:212)*10^3,FFTmag(2:end)/10^3)

legend('Rayleigh Integral Solution','k-Wave Solution')
hold off

%%

%% Plotting slice of xy plane %%
imagesc(x2*1e3,y2*1e3,PP(:,:,1)*1e-3);

shading flat
a = colorbar;
a.Label.String = 'Acoustic Pressure (kPa)';
axis equal
xlim([0 64])
ylim([-32 32])
title('Acoustic Pressure field from a bowl source in xy plane')
xlabel('Axis of propagation (mm)')
ylabel('Lateral position (mm)')


%%
%%% Error Analysis 
figure
vq2 = interp1(x2,pd,x_vec(1:276));
plot((x_vec(1:276)*10^3),vq2/10^3)
hold on 
plot((x_vec(1:276)*10^3),FFTmag(2:end)/10^3)
ylabel('Magnitude of Acoustic Pressure at Peak Frequency P(f) [kPa]')
xlabel('Axial Position (mm)')
title('Acoustic Pressure at Peak Frequency on Axial Positions after reaching steady-state with a bowl piston')
legend("Rayleigh's Integral",'k-Wave')
hold off
pctError = mean((abs(transpose(FFTmag(2:end))-vq2))./vq2*100)
D = (vq2-(transpose(FFTmag(2:end)))).^2;
S1= sum(D);
S2 = sum((vq2).^2);
ctest = 100*(sqrt(S1/S2))