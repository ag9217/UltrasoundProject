function sensor_data = kspaceFirstOrder3DG(varargin)
%KSPACEFIRSTORDER3DG 3D time-domain simulation of wave propagation on a GPU using C++ CUDA code.
%
% DESCRIPTION:
%     kspaceFirstOrder3DG provides an blind interface to the native C++
%     CUDA version of kspaceFirstOrder3D (called kspaceFirstOrder3D-CUDA)
%     in the same way as kspaceFirstOrder3DC.
%
%     This function requires the C++ binary/executable of 
%     kspaceFirstOrder3D-CUDA to be downloaded from
%     http://www.k-wave.org/download.php and placed in the "binaries"
%     directory of the k-Wave toolbox. 
% 
%     Note, not all input options are currently supported, and all
%     display options are ignored (only command line outputs are given).
%     See the k-Wave user manual for more information.
%
% USAGE:
%     see kspaceFirstOrder3DC
%
% ABOUT:
%     author          - Bradley Treeby
%     date            - 30th October 2013
%     last update     - 8th June 2017
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2013-2017 Bradley Treeby
%
% See also kspaceFirstOrder3D, kspaceFirstOrder3DC

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

% This function is essentially a wrapper and directly uses the capabilities
% of kspaceFirstOrder3DC by replacing the binary name with the name of the
% GPU binary. 

% check for a custom binary name
if any(strcmp('BinaryName', varargin))
    
    % if the binary name is given, directly pass to kspaceFirstOrder3DC
    sensor_data = kspaceFirstOrder3DC(varargin{:});
    
else

    % if the binary name is not given, specify to use the GPU binary
    if isunix
        binary_name = 'kspaceFirstOrder3D-CUDA';
    else
        binary_name = 'kspaceFirstOrder3D-CUDA.exe';
    end
    
    % pass this and the original inputs to kspaceFirstOrder3DC
    sensor_data = kspaceFirstOrder3DC(varargin{:}, 'BinaryName', binary_name);
    
end
