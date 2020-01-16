clearvars;

% simulation settings
DATA_CAST = 'single';

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 10;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 5;            % [grid points]

% set total number of grid points not including the PML
Nx = 128 - 2*PML_X_SIZE;    % [grid points]
Ny = 128 - 2*PML_Y_SIZE;    % [grid points]
Nz = 64 - 2*PML_Z_SIZE;     % [grid points]

% set desired grid size in the x-direction not including the PML
x = 40e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/Nx;                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
%% DEFINE THE MEDIUM PARAMETERS
% define the properties of the propagation medium
medium.sound_speed = 1540 * ones(Nx, Ny, Nz);      % [m/s]
medium.density = 1000 * ones(Nx, Ny, Nz);          % [kg/m^3]

% % 2nd denser layer
% medium.sound_speed(Nx/2:end,:,:) = 2000; 
% medium.density(Nx/2:end,:,:) = 5000;
% 
% medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
% medium.alpha_power = 1.5;
% medium.BonA = 6;


% create the time array
t_end = 40e-6;                  % [s]
kgrid.makeTime(medium.sound_speed, [], t_end);


%% DEFINE THE INPUT SIGNAL
% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 0.5e6;        % [Hz]
tone_burst_cycles = 10;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% =========================================================================
%% DEFINE THE ULTRASOUND TRANSDUCER
% physical properties of the transducer
transducer.number_elements = 22;    % total number of transducer elements
transducer.element_width = 1;       % width of each element [grid points]
transducer.element_length = 7;     % length of each element [grid points]
transducer.element_spacing = 0;     % spacing (kerf width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, 53, 29]);

% properties used to derive the beamforming delays
transducer.sound_speed = 1540;                  % sound speed [m/s]
transducer.focus_distance = 20e-3;              % focus distance [m]
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

%% Define cone
% Getting a logic array from the SolidWorks 3D model
[OUTPUTgrid] = VOXELISE(Nx/2,Ny/2,Nz,'Cone_model.stl','xyz');

OUTPUTgrid(Nx/2 + 1:Nx,:,:) = 0;

OUTPUTgrid(:,Ny/2 + 1:Ny,:) = 0;


% Type casting logic array to a different type 
Test = cast(OUTPUTgrid, 'double');

% Changing medium properties according to the shape of the model
% Also puts the cone on the centre of the Y-Z plane
for X = 1:length(Test(1:Nx,1,1))  
    for Y = 22:length(Test(1,1:Ny,1))
        for Z = 1:length(Test(1,1,1:Nz))
            if OUTPUTgrid(X,Y-21,Z) == 1
                 medium.density(X,Y,Z) = 10000;
                 medium.sound_speed(X,Y,Z) = 1000;
            end
        end
    end
end

% % Using function resize to change the resolution
% [Try] = cast(VOXELISE(64,128,128,'Cone_model.stl','xyz'),'single');
% mat_rs = resize(Try, [32 64 64]);
% 
% % Using the function VOXELISE but with lower resolution values for input
% [Low_res] = cast(VOXELISE(32,64,64,'Cone_model.stl','xyz'),'single');
% 
% % Initialize new 3D array that will join information from two different
% % resizing methods
% combined = cast(zeros(32, 64, 64), 'single');


% for X = 1:length(mat_rs(:,1,1))
%     for Y = 1:length(mat_rs(1,:,1))
%         for Z = 1:length(mat_rs(1,1,:))
%             if mat_rs(X,Y,Z) > 0 || Low_res(X,Y,Z) == 1
%                 combined(X,Y,Z) = 1;
%             end
% %             if mat_rs(X,Y,Z) > 0
% %                 mat_rs(X,Y,Z) = 1;
% %             end
%             
%         end
%     end
% end


% combined(33:128,:,:) = 0;
% combined(:,65:128,:) = 0;
% 
% for X = 1:length(combined(:,1,1))
%     for Y = 32:length(combined(1,:,1))
%         for Z = 1:length(combined(1,1,:))
%             if combined(X,Y-31,Z) == 1
%                 medium.density(X,Y,Z) = 10000;
%                 medium.sound_speed(X,Y,Z) = 1000;
%                 %Add other properties
%             end
%         end
%     end
% end

voxelPlot(medium.density);

%% DEFINE SENSOR MASK
% create a binary sensor mask with four detection positions
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask([Nx/4, Nx/2, 3*Nx/4], Ny/2, Nz/2) = 1;


%% RUN THE SIMULATION
% set the input settings
input_args = {'DisplayMask', transducer.all_elements_mask | sensor.mask, ...
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

