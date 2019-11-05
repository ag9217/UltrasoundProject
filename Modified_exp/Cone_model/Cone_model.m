clear all
Nx = 128;
Ny = 128;
Nz = 64;


[OUTPUTgrid] = VOXELISE(32,64,64,'Cone_model.stl','xyz');

% Filling the rest of the grid with zeros - has to match Nx,Ny,Nz
% dimensions
OUTPUTgrid(31:64,:,:) = 0;
OUTPUTgrid(:,65:128,:) = 0;
OUTPUTgrid(:,:,65:128) = 0;

Test = cast(OUTPUTgrid, 'uint8');

medium.sound_speed = 1000 * ones(Nx, Ny, Nz);	% [m/s]


for X = 1:length(Test(:,1,1))
    for Y = 1:length(Test(1,:,1))
        for Z = 1:length(Test(1,1,:))
            if OUTPUTgrid(X,Y,Z) == 1
                medium.sound_speed(X,Y,Z) = 2000;
            end
        end
    end
end


voxelPlot(medium.sound_speed)