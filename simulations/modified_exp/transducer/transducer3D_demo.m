% Group Project - 3rd Year                                     21/10/2019 
% Transducer part copied from bmode linear transducer example 
% 3D simulation coded around it to understand + play around with 
% specifics of transducer 
 
 

% create the computational grid
% set the size of the perfectly matched layer (PML)
% PML is a special type of anisotropic absorbing boundary layer
% by default this layer occupies a strip around the edge of the domain of
% 10 grid points in 3D
pml_x_size = 20;                % [grid points]
pml_y_size = 10;                % [grid points]
pml_z_size = 10;                % [grid points]
 

% set total number of grid points not including the PML
Nx = 256 - 2 * pml_x_size;      % [grid points]
Ny = 128 - 2 * pml_y_size;      % [grid points]
Nz = 128 - 2 * pml_z_size;      % [grid points]
dx = 50e-6;
dy = 50e-6;
dz = 50e-6;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
 

% define the medium properties
medium.sound_speed = 1500*ones(Nx, Ny, Nz);   % [m/s]
medium.sound_speed(1:50, :) = 1800;           % [m/s]
medium.density = 1040;                        % [kg/m^3]
kgrid.makeTime(medium.sound_speed);
 

% physical properties of the transducer
tr.number_elements = 5;     % total number of transducer elements
tr.element_width = 1;       % width of each element [grid points]
tr.element_length = 5;      % length of each element [grid points]
tr.element_spacing = 0;     % spacing (kerf  width) between the elements [grid points]
tr.radius = inf;            % radius of curvature of the transducer [m]
 

% calculate the width of the transducer in grid points
transducer_width = tr.number_elements * tr.element_width ...
    + (tr.number_elements - 1) * tr.element_spacing;
 

% use this to position the transducer in the middle of the computational grid
tr.position = round([1, Ny/2 - transducer_width/2, Nz/2 - tr.element_length/2]);
 

% properties used to derive the beamforming delays
%tr.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
tr.sound_speed = 1540;                   % sound speed [m/s]
tr.focus_distance = 20e-3;               % focus distance [m]
tr.steering_angle = 10;                  % steering angle [degrees]
 

% apodization (optical filter technique)
tr.transmit_apodization = 'Hanning';    
tr.receive_apodization = 'Rectangular';
 

% define the transducer elements that are currently active
number_active_elements = 10;
tr.active_elements = ones(tr.number_elements, 1);

source_strength = 1e6;
tone_burst_freq = 1.5e6;
tone_burst_cycles = 4;
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% append input signal used to drive the transducer
tr.input_signal = input_signal;
 

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, tr);
 

% run a simulation using the same transducer as both source and sensor
sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer);
    

