clear all;

% % Grid dimensions 
% Nx = 128;
% Ny = 128;
% Nz = 64;

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 10;            % [grid points]

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

% Getting a logic array from the SolidWorks 3D model
[OUTPUTgrid] = VOXELISE(32,64,64,'Cone_model.stl','xyz');


% Filling the rest of the grid with zeros - has to match Nx,Ny,Nz
% dimensions
OUTPUTgrid(31:128,:,:) = 0;
OUTPUTgrid(:,65:128,:) = 0;

% Type casting logic array to a different type 
Test = cast(OUTPUTgrid, 'double');

% Initializing uniform mediums speed throught entire grid
medium.sound_speed = 1000 * ones(Nx, Ny, Nz);	% [m/s]

% Changing medium properties according to the shape of the model
% Also puts the cone on the centre of the Y-Z plane
for X = 1:length(Test(:,1,1))  
    for Y = 32:length(Test(1,:,1))
        for Z = 1:length(Test(1,1,:))
            if OUTPUTgrid(X,Y-31,Z) == 1
                medium.sound_speed(X,Y,Z) = 2000;
            end
        end
    end
end

voxelPlot(medium.sound_speed)
title('Cone in low res - unmodified')
%% Trying to find better voxelplot when lower resolution is used
%                           WORKS!! ?
% Idea have to different ways of reducing the resolution of the cone
% Join them together so that the loss of information is minimal

% Using function resize to change the resolution
[Try] = cast(VOXELISE(64,128,128,'Cone_model.stl','xyz'),'single');
mat_rs = resize(Try, [32 64 64]);

% Using the function VOXELISE but with lower resolution values for input
[Low_res] = cast(VOXELISE(32,64,64,'Cone_model.stl','xyz'),'single');

% Initialize new 3D array that will join information from two different
% resizing methods
combined = cast(zeros(32, 64, 64), 'single');


for X = 1:length(mat_rs(:,1,1))
    for Y = 1:length(mat_rs(1,:,1))
        for Z = 1:length(mat_rs(1,1,:))
            if mat_rs(X,Y,Z) > 0 || Low_res(X,Y,Z) == 1
                combined(X,Y,Z) = 1;
            end
%             if mat_rs(X,Y,Z) > 0
%                 mat_rs(X,Y,Z) = 1;
%             end
            
        end
    end
end


combined(33:128,:,:) = 0;
combined(:,65:128,:) = 0;

medium.density = 1000 * ones(Nx, Ny, Nz);
for X = 1:length(combined(:,1,1))
    for Y = 32:length(combined(1,:,1))
        for Z = 1:length(combined(1,1,:))
            if combined(X,Y-31,Z) == 1
                medium.density(X,Y,Z) = 2000;
            end
        end
    end
end

voxelPlot(medium.density)
title('Cone in low res - modified')

