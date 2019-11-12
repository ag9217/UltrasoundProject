clear all;

% Grid dimensions 
Nx = 128;
Ny = 128;
Nz = 64;

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
            if mat_rs(X,Y,Z) > 0
                mat_rs(X,Y,Z) = 1;
            end
            
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

