%% Combining transducer and cone
clear all;
close all;
% Simulation settings
DATA_CAST = 'single';

% Defining the kgrid
% PML
PML_X_SIZE = 10;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 20;            % [grid points]

% Absolute grid values
Nx = 64;    % [grid points]
Ny = 128;    % [grid points]
Nz = 128;     % [grid points]

% Defining size of each element in the kgrid
dx = 2e-3; 
dy = 2e-3;
dz = 2e-3;

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

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
[Cone] = VOXELISE(cone_gridPoints_Z,cone_gridPoints_X,cone_gridPoints_Y,'Cone_model.stl','xyz');
Cone = cast(Cone, 'single');

% Match cone dimensions with kgrid - cone located at the corner
% Warning: coordinets can be a bit confusing
Cone(cone_gridPoints_Z+1:Nx,:,:) = 0;
Cone(:,cone_gridPoints_X+1:Nz,:) = 0;
Cone(:,:,cone_gridPoints_Y+1:Ny) = 0;


% Shifting the cone to the middle of the grid
difference = ((total_length_Y - R)/2)/(dx*1000);
Cone = circshift(Cone, [0 difference difference]);

% Swapping the dimensions of the cone to match kgrid
%Cone = permute(Cone, [2 3 1]); % second dimension becomes first, third one becomes second and first one becomes third

% Defining medium properties - look up real values
medium.sound_speed = 1540 * ones(Nx,Ny,Nz);
medium.density = 1000 * ones (Nx,Ny,Nz);

medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% Changing sounds properties according to the model
for X = 1:length(Cone(:,1,1))  
    for Y = 1:length(Cone(1,:,1))
        for Z = 1:length(Cone(1,1,:))
            if Cone(X,Y,Z) == 1
                medium.sound_speed(X,Y,Z) = 2000;
                medium.density(X,Y,Z) = 5000;
            end
        end
    end
end


%% Defining Transducer

% create the time array
t_end = 40e-6;                  % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);


% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 0.5e6;        % [Hz]
tone_burst_cycles = 1;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% physical properties of the transducer
transducer.number_elements = 72;    % total number of transducer elements
transducer.element_width = 1;       % width of each element [grid points]
transducer.element_length = 12;     % length of each element [grid points]
transducer.element_spacing = 0;     % spacing (kerf width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
%Positioning the transducer on the middle point of the Y-Z plane
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);


% properties used to derive the beamforming delays
transducer.sound_speed = 1540;                  % sound speed [m/s]
transducer.focus_distance = 40e-3;              % focus distance [m]
transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = zeros(transducer.number_elements, 1);
transducer.active_elements(21:52) = 1;

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;


% create voxel plot of transducer mask and cone
voxelPlot(single(transducer.active_elements_mask | Cone));

%% Sensor mask
% create a binary sensor mask with four detection positions
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask([Nx/4, Nx/2, 3*Nx/4], Ny/2, Nz/2) = 1;

%% SIMULATION

input_args = {'DisplayMask', transducer.all_elements_mask | sensor.mask ...
    'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    'DataCast', DATA_CAST, 'PlotScale', [-1/2, 1/2] * source_strength};

% run the simulation
[sensor_data] = kspaceFirstOrder3D(kgrid, medium, transducer, sensor, input_args{:});


% calculate the amplitude spectrum of the input signal and the signal
% recorded each of the sensor positions 
[f_input, as_input] = spect([input_signal, zeros(1, 2 * length(input_signal))], 1/kgrid.dt);
[~, as_1] = spect(sensor_data(1, :), 1/kgrid.dt);
[~, as_2] = spect(sensor_data(2, :), 1/kgrid.dt);
[f, as_3] = spect(sensor_data(3, :), 1/kgrid.dt);