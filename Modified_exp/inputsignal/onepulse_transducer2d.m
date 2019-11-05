% Try to make the source (transducer) a single pulse             05/11/2019
% instead of a continuous pulse, want to achieve this to see 
% if we can perceive the reflected wave
% ATTEMPT 2 (using kspaceFirstOrder --> transducer)

close all
clear all

% create computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.1e-3;        % grid point spacing in the x direction [m]
dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define properties medium
medium.sound_speed              = 1500 * ones(Nx, Ny);
medium.density                  = 1000 * ones(Nx, Ny);   % [kg/m^3]

% fibroid
fibroid = makeDisc(Nx, Ny, 100, 64, 20, true);
medium.sound_speed(fibroid==1)=1900;
medium.density(fibroid==1) = 10000;

kgrid.makeTime(medium.sound_speed);

% create initial pressure distribution using makeDisc
%disc_magnitude = 5; % [Pa]
%disc_x_pos = 30;    % [grid points]
%disc_y_pos = 64;    % [grid points]
%disc_radius = 5;    % [grid points]
%source.p0 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

% define a curved transducer element
arc_pos = [20, 20];         % [grid points]    
radius = 60;                % [grid points]
diameter = 81;              % [grid points]
focus_pos = [Nx/2, Nx/2];   % [grid points]
source.p_mask = makeArc([Nx, Ny], arc_pos, radius, diameter, focus_pos);

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 1.5e6;        % [Hz]
tone_burst_cycles = 2;
%tone_burst_cycles = 1;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
source.p = input_signal;

% filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);
% create a display mask to display the transducer
display_mask = source.p_mask;


% define sensor
%sensor.mask = [1, 1, Nx, Ny].';
sensor.mask = makeCircle(Nx, Ny, Nx/2, Ny/2, 20);

% run simulation 
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% reorder the simulation data
sensor_data_reordered = reorderSensorData(kgrid, sensor, sensor_data);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot source, sensor, and position of the interface
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, ...
    double(source.p_mask | sensor.mask), [-1 1]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% plot the re-ordered sensor data
figure;
imagesc(sensor_data_reordered, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;




