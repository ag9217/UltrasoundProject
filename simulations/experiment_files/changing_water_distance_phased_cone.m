% Simulating B-mode Ultrasound Images Example
%
% This example illustrates how k-Wave can be used for the simulation of
% B-mode ultrasound images using a phased-array or sector transducer. It
% builds on the Simulating B-mode Ultrasound Images Example.
%
% To allow the simulated scan line data to be processed multiple times with
% different settings, the simulated RF data is saved to disk. This can be
% reloaded by setting RUN_SIMULATION = false within the example m-file. The
% data can also be downloaded from
% http://www.k-wave.org/datasets/example_us_phased_array_scan_lines.mat 
%
% author: Bradley Treeby
% date: 7th September 2012
% last update: 6th August 2017
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

%#ok<*UNRCH>

clearvars;

% simulation settings
DATA_CAST       = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 15;                % [grid points]
pml_y_size = 10;                % [grid points]
pml_z_size = 10;                % [grid points]

% set total number of grid points not including the PML
sc = 1;
Nx = 256/sc - 2*pml_x_size;     % [grid points]
Ny = 128/sc - 2*pml_y_size;     % [grid points]
Nz = 128/sc - 2*pml_z_size;     % [grid points]

% set desired grid size in the x-direction not including the PML
x = 150e-3;                     % [m]

% calculate the spacing between the grid points
dx = x / Nx;                    % [m]
dy = dx;                        % [m]
dz = dx;                        % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% changing dimensions to allow lower resolution cone

% set total number of grid points not including the PML
Nx = 128;                       % [grid points]
Ny = 128;                       % [grid points]
Nz = 128;                       % [grid points]

% set desired grid size in the x-direction not including the PML
x = 40e-3;                      % [m]

% calculate the spacing between the grid points
dx = x / Nx;                    % [m]
dy = dx;                        % [m]
dz = dx;                        % [m]

% Cone dimensions
R = 220; % Outer diameter [mm]
t = 50; % Thickness of the cone [mm] - check with Solidworks

% Calculating total length of the kgrid in mm
total_length_X = Nx*dx*1000; % [mm]
total_length_Y = Ny*dy*1000; % [mm]
total_length_Z = Nz*dz*1000; % [mm]

% Size of the grid that cone is occupying

cone_gridPoints_X = round(R/(dx*1000));
cone_gridPoints_Y = round(R/(dy*1000));
cone_gridPoints_Z = round(t/(dz*1000));

% Importing the cone
[Cone] = VOXELISE(cone_gridPoints_Z/16,cone_gridPoints_X/16,cone_gridPoints_Y/16,'Cone_model.stl','xyz');
Cone = cast(Cone, 'single');
% Match cone dimensions with kgrid - cone located at the corner
% Warning: coordinets can be a bit confusing
Cone(cone_gridPoints_Z+1:Nx,:,:) = 0;
Cone(:,cone_gridPoints_X+1:Nz,:) = 0;
Cone(:,:,cone_gridPoints_Y+1:Ny) = 0;
Cone = padarray(Cone, [0,Ny/4,Nz/4], 'both');
Cone = circshift(Cone, [0 0 -1]);
voxelPlot(Cone)
Cone(1,:,:) = 0;
Cone(2,:,:) = 0;
% Shifting the cone to the middle of the grid

% setting dimensions back to as specified above

% set total number of grid points not including the PML
sc = 1;
Nx = 256/sc - 2*pml_x_size;     % [grid points]
Ny = 128/sc - 2*pml_y_size;     % [grid points]
Nz = 128/sc - 2*pml_z_size;     % [grid points]

% set desired grid size in the x-direction not including the PML
x = 150e-3;                     % [m]

% calculate the spacing between the grid points
dx = x / Nx;                    % [m]
dy = dx;                        % [m]
dz = dx;                        % [m]

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
medium.alpha_coeff = 0.1;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 0.5e6 / sc;     % [Hz]
tone_burst_cycles = 4;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (c0 * rho0)) .* input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% define the physical properties of the phased array transducer
transducer.number_elements = 16 / sc;       % total number of transducer elements
transducer.element_width = 1;               % width of each element [grid points]
transducer.element_length = 20 / sc;        % length of each element [grid points]
transducer.element_spacing = 0;             % spacing (kerf  width) between the elements [grid points]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([7, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = c0;                    % sound speed [m/s]
transducer.focus_distance = 30e-3;              % focus distance [m]
transducer.elevation_focus_distance = 30e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]
transducer.steering_angle_max = 32;             % maximum steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Hanning';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = ones(transducer.number_elements, 1);

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% =========================================================================
% DEFINE THE MEDIUM PROPERTIES
% =========================================================================

% range of steering angles to test
steering_angles = -2:2:2;

% preallocate the storage
number_scan_lines = length(steering_angles);
scan_lines = zeros(number_scan_lines, kgrid.Nt);

% define a random distribution of scatterers for the medium
background_map_mean = 1;
background_map_std = 0.008;
background_map = background_map_mean + background_map_std * randn([Nx, Ny, Nz]);

% define a random distribution of scatterers for the highly scattering
% region
scattering_map = randn([Nx, Ny, Nz]);
scattering_c0 = c0 + 25 + 75 * scattering_map;
scattering_c0(scattering_c0 > 1600) = 1600;
scattering_c0(scattering_c0 < 1400) = 1400;
scattering_rho0 = scattering_c0 / 1.5;

% define properties
sound_speed_map = c0 * ones(Nx, Ny, Nz) .* background_map;
density_map = rho0 * ones(Nx, Ny, Nz) .* background_map;

% defining water layer 20 degrees C
water_dist = 10e-3;
water_layer = round(water_dist/dx);
sound_speed_map(1:water_layer,:,:) = 1481;
density_map(1:water_layer,:,:) = 998;
% reapplying randomness to newly defined layer
sound_speed_map(1:water_layer,:,:) = sound_speed_map(1:water_layer,:,:) .* background_map(1:water_layer,:,:);

% defining skin layer (1mm)
skin_layer = 1e-3/dx;
water_skin = water_layer + skin_layer;
sound_speed_map(water_layer:water_skin,:,:) = 1624;
density_map(water_layer:water_skin,:,:) = 1050;
sound_speed_map(water_layer:water_skin,:,:) = sound_speed_map(water_layer:water_skin,:,:) .* background_map(water_layer:water_skin,:,:);

% defining subcut fat (20mm)
fat_layer = 20e-3/dx;
water_skin_fat = water_skin + fat_layer;
sound_speed_map(water_skin:water_skin_fat,:,:) = 1450;
density_map(water_skin:water_skin_fat,:,:) = 900;
sound_speed_map(water_skin:water_skin_fat,:,:) = sound_speed_map(water_skin:water_skin_fat,:,:) .* background_map(water_skin:water_skin_fat,:,:);

% defining muscle (10mm)
muscle_layer = 10e-3/dx;
water_skin_fat_musc = water_skin_fat + muscle_layer;
sound_speed_map(water_skin_fat:water_skin_fat_musc,:,:) = 1580;
density_map(water_skin_fat:water_skin_fat_musc,:,:) = 1100;
sound_speed_map(water_skin_fat:water_skin_fat_musc,:,:) = sound_speed_map(water_skin_fat:water_skin_fat_musc,:,:) .* background_map(water_skin_fat:water_skin_fat_musc,:,:);

% define a sphere for a highly scattering region
Nx_tot = Nx;
Ny_tot = Ny + number_scan_lines * transducer.element_width;
Nz_tot = Nz;

radius = 5e-3/dx;
x_pos = water_layer + skin_layer + fat_layer + muscle_layer + radius;
y_pos = Ny_tot/2;
z_pos = Nz_tot/2;
% water_layer + round((fib_tiss_distance + radius)/dx)
scattering_region2 = makeBall(Nx_tot, Ny, Nz_tot, x_pos, y_pos, z_pos, radius);

% assign region
sound_speed_map(scattering_region2 == 1) = scattering_c0(scattering_region2 == 1);
density_map(scattering_region2 == 1) = scattering_rho0(scattering_region2 == 1);

% assign to the medium inputs
medium.sound_speed = sound_speed_map;
medium.density = density_map;

%cone material properties
for X = 1:length(Cone(:,1,1))  
    for Y = 1:length(Cone(1,:,1))
        for Z = 1:length(Cone(1,1,:))
            if Cone(X,Y,Z) == 1
                sound_speed_map(X,Y,Z) = 5800;
                density_map(X,Y,Z) = 7700;
            end
        end
    end
end

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

% run the simulation if set to true, otherwise, load previous results
if RUN_SIMULATION
    
    % loop through the range of angles to test
    for angle_index = 1:number_scan_lines
        
        % update the command line status
        disp('');
        disp(['Computing scan line ' num2str(angle_index) ' of ' num2str(number_scan_lines)]);

        % update the current steering angle
        transducer.steering_angle = steering_angles(angle_index);
        
        % run the simulation
        sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer, transducer, input_args{:});

        % extract the scan line from the sensor data
        scan_lines(angle_index, :) = transducer.scan_line(sensor_data);
        
    end

    % save the scan lines to disk
    save example_us_phased_array_scan_lines scan_lines;
    
else
    
    % load the scan lines from disk
    load example_us_phased_array_scan_lines
    
end

% trim the delay offset from the scan line data
t0_offset = round(length(input_signal) / 2) + (transducer.appended_zeros - transducer.beamforming_delays_offset);
scan_lines = scan_lines(:, t0_offset:end);

% get the new length of the scan lines
Nt = length(scan_lines(1, :));

% =========================================================================
% PROCESS THE RESULTS
% =========================================================================

% -----------------------------
% Remove Input Signal
% -----------------------------

% create a window to set the first part of each scan line to zero to remove
% interference from the input signal
scan_line_win = getWin(Nt * 2, 'Tukey', 'Param', 0.05).';
scan_line_win = [zeros(1, t0_offset * 2), scan_line_win(1:end/2 - t0_offset * 2)];

% apply the window to each of the scan lines
scan_lines = bsxfun(@times, scan_line_win, scan_lines);

% -----------------------------
% Time Gain Compensation
% -----------------------------

% create radius variable
r = c0 * (1:Nt) * kgrid.dt / 2;    % [m]

% create time gain compensation function based on attenuation value,
% transmit frequency, and round trip distance
tgc_alpha = 0.4;       % [dB/(MHz cm)]
tgc = exp(2 * tgc_alpha * tone_burst_freq * 1e-6 * r * 100);

% apply the time gain compensation to each of the scan lines
scan_lines = bsxfun(@times, tgc, scan_lines);

% -----------------------------
% Frequency Filtering
% -----------------------------

% filter the scan lines using both the transmit frequency and the second
% harmonic
scan_lines_fund = gaussianFilter(scan_lines, 1/kgrid.dt, tone_burst_freq, 100, true);
scan_lines_harm = gaussianFilter(scan_lines, 1/kgrid.dt, 2 * tone_burst_freq, 30, true);

% -----------------------------
% Envelope Detection
% -----------------------------

% envelope detection
scan_lines_fund = envelopeDetection(scan_lines_fund);
scan_lines_harm = envelopeDetection(scan_lines_harm);

% -----------------------------
% Log Compression
% -----------------------------

% normalised log compression
compression_ratio = 3;
scan_lines_fund = logCompression(scan_lines_fund, compression_ratio, true);
scan_lines_harm = logCompression(scan_lines_harm, compression_ratio, true);

% -----------------------------
% Scan Conversion
% -----------------------------

% set the desired size of the image
image_size = [Nx * dx, Ny * dy];

% convert the data from polar coordinates to Cartesian coordinates for
% display
b_mode_fund = scanConversion(scan_lines_fund, steering_angles, image_size, c0, kgrid.dt);
b_mode_harm = scanConversion(scan_lines_harm, steering_angles, image_size, c0, kgrid.dt);

% =========================================================================
% VISUALISATION
% =========================================================================

% create the axis variables
x_axis = [0, Nx * dx * 1e3];    % [mm]
y_axis = [0, Ny * dy * 1e3];    % [mm]

% plot the data before and after scan conversion
figure;
subplot(1, 3, 1);
imagesc(steering_angles, x_axis, scan_lines.');
axis square;
xlabel('Steering angle [deg]');
ylabel('Depth [mm]');
title('Raw Scan-Line Data');

subplot(1, 3, 2);
imagesc(steering_angles, x_axis, scan_lines_fund.');
axis square;
xlabel('Steering angle [deg]');
ylabel('Depth [mm]');
title('Processed Scan-Line Data');

subplot(1, 3, 3);
imagesc(y_axis, x_axis, b_mode_fund);
axis square;
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');
title('B-Mode Image');
colormap(gray);

scaleFig(2, 1);

% plot the medium and the B-mode images
figure;
subplot(1, 3, 1);
imagesc(y_axis, x_axis, medium.sound_speed(:, :, end/2));
axis image;
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');
title('Scattering Phantom');

subplot(1, 3, 2);
imagesc(y_axis, x_axis, b_mode_fund);
axis image;
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');
title('B-Mode Image');

subplot(1, 3, 3);
imagesc(y_axis, x_axis, b_mode_harm);
colormap(gray);
axis image;
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');
title('Harmonic Image');

scaleFig(2, 1);

% =========================================================================
% VISUALISATION OF SIMULATION LAYOUT
% =========================================================================

% physical properties of the transducer
transducer_plot.number_elements = 16 / sc;       % total number of transducer elements
transducer_plot.element_width = 1;               % width of each element [grid points]
transducer_plot.element_length = 20 / sc;        % length of each element [grid points]
transducer_plot.element_spacing = 0;             % spacing (kerf  width) between the elements [grid points]
transducer_plot.radius = inf;

% transducer position
transducer_plot.position = round([7, Ny_tot/2 - transducer_width/2, Nz_tot/2 - transducer.element_length/2]);

% create expanded grid
kgrid_plot = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
kgrid_plot.setTime(1, 1);

% create the transducer using the defined settings
transducer_plot = kWaveTransducer(kgrid_plot, transducer_plot);
hold on;
% create voxel plot of transducer mask and
out = sound_speed_map > 5000;
voxelPlot(single(transducer_plot.active_elements_mask | scattering_region2 | out));
view(26, 48);