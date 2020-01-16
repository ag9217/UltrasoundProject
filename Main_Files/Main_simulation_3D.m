% NOTE: This code will NOT work on your computer unless you have a CUDA
% enabled graphics card. If you wish to run the simulation without a CUDA
% accelerated GPU (I would not recommend it), set DATA_CAST = 'single' and
% remove the G from kspaceFirstOrder3DG (line 179)

clearvars;
close all;

% simulation settings
DATA_CAST       = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = true;         % set to false to reload previous results instead of running simulation

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 20;                % [grid points]
pml_y_size = 10;                % [grid points]
pml_z_size = 10;                % [grid points]

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

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% CONE DIMENSIONS
% =========================================================================

% x = 40e-3;                      % [m]
% Nx = 128;
% dx = x / Nx;                    % [m]
% dy = dx;                        % [m]
% dz = dx; 

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
[Cone] = VOXELISE(cone_gridPoints_Z/8,cone_gridPoints_X/8,cone_gridPoints_Y/8,'Cone_model.stl','xyz');
Cone = cast(Cone, 'single');
% Match cone dimensions with kgrid - cone located at the corner
% Warning: coordinets can be a bit confusing
Cone(cone_gridPoints_Z+1:Nx,:,:) = 0;
Cone(:,cone_gridPoints_X+1:Nz,:) = 0;
Cone(:,:,cone_gridPoints_Y+1:Ny) = 0;
Cone = padarray(Cone, [0,(160-88)/2,(128-88)/2], 'both');
Cone = circshift(Cone, [0 0 -1]);
voxelPlot(Cone)
Cone(1,:,:) = 0;
Cone(2,:,:) = 0;
% Shifting the cone to the middle of the grid

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% x = 100e-3;                      % [m]
% Nx = 256;
% dx = x / Nx;                    % [m]
% dy = dx;                        % [m]
% dz = dx; 

% define the properties of the propagation medium
% these properties are very similair to human tissue (fat)
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, 0.01, t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 1.5e6;        % [Hz]
tone_burst_cycles = 4;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (c0 * rho0)) .* input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 8;  	% total number of transducer elements
transducer.element_width = 0.5;       % width of each element [grid points]
transducer.element_length = 24;  	% length of each element [grid points]
transducer.element_spacing = 0;  	% spacing (kerf  width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = c0;                    % sound speed [m/s]
transducer.focus_distance = 20e-3;              % focus distance [m]
transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Hanning';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
number_active_elements = 8;
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

% define a large image size to move across
number_scan_lines = 32;
Nx_tot = Nx;
Ny_tot = Ny + number_scan_lines * transducer.element_width;
Nz_tot = Nz;

% define a random distribution of scatterers for the medium
background_map_mean = 1;
background_map_std = 0.008;
% uniform distribution of medium scattering
background_map = background_map_mean + background_map_std * randn([Nx_tot, Ny_tot, Nz_tot]);

% define a random distribution of scatterers for the highly scattering
% region
scattering_map = randn([Nx_tot, Ny_tot, Nz_tot]);
scattering_c0 = c0 + 25 + 75 * scattering_map;
scattering_c0(scattering_c0 > 1600) = 1600;
scattering_c0(scattering_c0 < 1400) = 1400;
scattering_rho0 = scattering_c0 / 1.5;

% define properties
sound_speed_map = c0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;
density_map = rho0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;

% ###### defining water layer ######
sound_speed_map(1:Nx/2,:,:) = 1500;
% reapplying randomness to newly defined layer
sound_speed_map(1:Nx/2,:,:) = sound_speed_map(1:Nx/2,:,:) .* background_map(1:Nx/2,:,:);

% define a sphere for a highly scattering region
radius = 6e-3;      % [m]
x_pos = 27.5e-3;    % [m]
y_pos = 20.5e-3;      % [m]
scattering_region2 = makeBall(Nx_tot, Ny_tot, Nz_tot, Nx_tot - round(radius/(x/128)), Ny_tot/2, Nz_tot/2, round(radius/(x/64)));

%cone material properties
% Importing model into the kgrid
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

% assign region
sound_speed_map(scattering_region2 == 1) = scattering_c0(scattering_region2 == 1);
density_map(scattering_region2 == 1) = scattering_rho0(scattering_region2 == 1);

%sound_speed_map(1:Nx/2, :,:) = 1480;

% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% preallocate the storage
scan_lines = zeros(number_scan_lines, kgrid.Nt);

% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

% run the simulation if set to true, otherwise, load previous results from
% disk
if RUN_SIMULATION

    % set medium position
    medium_position = 1;
    
    % loop through the scan lines
    for scan_line_index = 1:number_scan_lines
        
        % update the command line status
        disp('');
        disp(['Computing scan line ' num2str(scan_line_index) ' of ' num2str(number_scan_lines)]);

        % load the current section of the medium
        medium.sound_speed = sound_speed_map(:, medium_position:medium_position + Ny - 1, :);
        medium.density = density_map(:, medium_position:medium_position + Ny - 1, :);
        
        % run the simulation
        sensor_data = kspaceFirstOrder3DG(kgrid, medium, transducer, transducer, input_args{:});

        % extract the scan line from the sensor data
        scan_lines(scan_line_index, :) = transducer.scan_line(sensor_data);
          
        % update medium position
        medium_position = medium_position + transducer.element_width;

        end

    % save the scan lines to disk
    save example_us_bmode_scan_lines scan_lines;
    
    else
    
    % load the scan lines from disk
    load example_us_bmode_scan_lines;
    
    end

% =========================================================================
% PROCESS THE RESULTS
% =========================================================================
% Results are processed as these basic steps are also performed by a modern
% diagnostic ultrasound scanner
% -----------------------------
% Remove Input Signal
% -----------------------------

% create a window to set the first part of each scan line to zero to remove
% interference from the input signal
scan_line_win = getWin(kgrid.Nt * 2, 'Tukey', 'Param', 0.05).';
scan_line_win = [zeros(1, length(input_signal) * 2), scan_line_win(1:end/2 - length(input_signal) * 2)];

% apply the window to each of the scan lines
scan_lines = bsxfun(@times, scan_line_win, scan_lines);

% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(1, :) = scan_lines(end/2, :);

% -----------------------------
% Time Gain Compensation
% -----------------------------
% Reflections from deeper tissue features will appear weaker 
% create radius variable assuming that t0 corresponds to the middle of the
% input signal
t0 = length(input_signal) * kgrid.dt / 2;
r = c0 * ( (1:length(kgrid.t_array)) * kgrid.dt / 2 - t0);    % [m]

% create time gain compensation function based on attenuation value,
% transmit frequency, and round trip distance
tgc_alpha = 0.4;       % [dB/(MHz cm)]
tgc = exp(2 * tgc_alpha * tone_burst_freq * 1e-6 * r * 100);

% apply the time gain compensation to each of the scan lines
scan_lines = bsxfun(@times, tgc, scan_lines);

% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(2, :) = scan_lines(end/2, :);

% -----------------------------
% Frequency Filtering
% -----------------------------

% filter the scan lines using both the transmit frequency and the second
% harmonic
scan_lines_fund = gaussianFilter(scan_lines, 1/kgrid.dt, tone_burst_freq, 100, true);
set(gca, 'XLim', [0, 6]);
scan_lines_harm = gaussianFilter(scan_lines, 1/kgrid.dt, 2 * tone_burst_freq, 30, true);
set(gca, 'XLim', [0, 6]);

% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(3, :) = scan_lines_fund(end/2, :);

% -----------------------------
% Envelope Detection
% -----------------------------

% envelope detection
scan_lines_fund = envelopeDetection(scan_lines_fund);
scan_lines_harm = envelopeDetection(scan_lines_harm);

% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(4, :) = scan_lines_fund(end/2, :);

% -----------------------------
% Log Compression
% -----------------------------

% normalised log compression
compression_ratio = 3;
scan_lines_fund = logCompression(scan_lines_fund, compression_ratio, true);
scan_lines_harm = logCompression(scan_lines_harm, compression_ratio, true);

% store a copy of the middle scan line to illustrate the effects of each
% processing step
scan_line_example(5, :) = scan_lines_fund(end/2, :);

% -----------------------------
% Scan Conversion
% -----------------------------

% upsample the image using linear interpolation
scale_factor = 2;
scan_lines_fund = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_fund, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');
scan_lines_harm = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_harm, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the medium, truncated to the field of view
figure;
imagesc((0:number_scan_lines * transducer.element_width - 1) * dy * 1e3, (0:Nx_tot-1) * dx * 1e3, sound_speed_map(:, 1 + Ny/2:end - Ny/2, Nz/2));
axis image;
colormap(gray);
set(gca, 'YLim', [5, 40]);
title('Scattering Phantom');
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');

% plot the processing steps
figure;
%stackedPlot(kgrid.t_array * 1e6, {'1. Beamformed Signal', '2. Time Gain Compensation', '3. Frequency Filtering', '4. Envelope Detection', '5. Log Compression'}, scan_line_example);
xlabel('Time [\mus]');
set(gca, 'XLim', [5, t_end * 1e6]);

% plot the processed b-mode ultrasound image
figure;
horz_axis = (0:length(scan_lines_fund(:, 1)) - 1) * transducer.element_width * dy / scale_factor * 1e3;
imagesc(horz_axis, r * 1e3, scan_lines_fund.');
axis image;
colormap(gray);
set(gca, 'YLim', [5, 40]);
title('B-mode Image');
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');

% plot the processed harmonic ultrasound image
figure;
imagesc(horz_axis, r * 1e3, scan_lines_harm.');
axis image;
colormap(gray);
set(gca, 'YLim', [5, 40]);
title('Harmonic Image');
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');

% =========================================================================
% VISUALISATION OF SIMULATION LAYOUT
% =========================================================================

% uncomment to generate a voxel plot of the simulation layout

% physical properties of the transducer
transducer_plot.number_elements = 8 + number_scan_lines - 1;
transducer_plot.element_width = 0.5;
transducer_plot.element_length = 24;
transducer_plot.element_spacing = 0;
transducer_plot.radius = inf;

% transducer position
transducer_plot.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% create expanded grid
kgrid_plot = kWaveGrid(Nx_tot, dx, Ny_tot, dy, Nz, dz);
kgrid_plot.setTime(1, 1);

% create the transducer using the defined settings
transducer_plot = kWaveTransducer(kgrid_plot, transducer_plot);
hold on;
% create voxel plot of transducer mask and
out = sound_speed_map > 5000;
voxelPlot(single(transducer_plot.active_elements_mask | scattering_region2 | out));
view(26, 48);
