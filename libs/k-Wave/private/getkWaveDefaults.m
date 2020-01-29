% DESCRIPTION:
%     Subscript to set the literals and defaults used in the fluid and
%     elastic simulation codes.
%
% ABOUT:
%     author      - Bradley Treeby
%     date        - 11th February 2014
%     last update - 8th June 2017
%       
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2014-2017 Bradley Treeby

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

%#ok<*NASGU>

% =========================================================================
% SIMULATION FUNCTIONS
% =========================================================================

% number of input variables required to run the simulation codes
NUM_REQ_INPUT_VARIABLES             = 4;  

% set defaults for literals that can be changed using optional inputs
% (these are used in all codes)
CARTESIAN_INTERP_DEF                = 'linear';
CREATE_LOG_DEF                      = false;
DATA_CAST_DEF                       = 'off';
DATA_RECAST_DEF                     = false;
DISPLAY_MASK_DEF                    = 'default';
LOG_SCALE_DEF                       = false;
LOG_SCALE_COMPRESSION_FACTOR_DEF    = 0.02;
MOVIE_ARGS_DEF                      = {};
PLOT_FREQ_DEF                       = 10;
PLOT_LAYOUT_DEF                     = false;
PLOT_SIM_DEF                        = true;
PLOT_PML_DEF                        = true;
PML_ALPHA_DEF                       = 2;
PML_INSIDE_DEF                      = true;
RECORD_MOVIE_DEF                    = false;
SCALE_SOURCE_TERMS_DEF              = true;
SMOOTH_P0_DEF                       = true;
SMOOTH_C0_DEF                       = false;
SMOOTH_RHO0_DEF                     = false;
SOURCE_S_MODE_DEF                   = 'additive';
SOURCE_P_MODE_DEF                   = 'additive';
SOURCE_U_MODE_DEF                   = 'additive';
USE_KSPACE_DEF                      = true;
USE_SG_DEF                          = true;

% set defaults for literals that depend on predefined values
if exist('MFILE', 'var')
    MOVIE_NAME_DEF                  = [getDateString, '-', MFILE];
else
    MOVIE_NAME_DEF                  = getDateString;
end

% set defaults for literals that can be changed using optional inputs
% (these are only used in 1D)
USE_FINITE_DIFFERENCE_DEF           = false;

% set defaults for literals that can be changed using optional inputs
% (these are only used in 2D)
MESH_PLOT_DEF                       = false;
FORCE_TSEARCH                       = false;
DIRECTIVITY_PATTERN_DEF             = 'pressure';
DIRECTIVITY_SIZE_SCALE_FACTOR_DEF   = 10;

% set defaults for literals that can be changed using optional inputs
% (these are only used in 3D)
SAVE_TO_DISK_DEF                    = false;
SAVE_TO_DISK_FILENAME_DEF           = 'kwave_input_data.h5';
SAVE_TO_DISK_EXIT_DEF               = true;
STREAM_TO_DISK_DEF                  = false;
STREAM_TO_DISK_STEPS_DEF            = 200;
STREAM_TO_DISK_FILENAME             = 'temp_sensor_data.bin';

% set default movie compression
MOVIE_PROFILE_DEF                   = 'Uncompressed AVI';

% set additional literals that can't be changed using optional inputs
COLOR_MAP                           = getColorMap;
ESTIMATE_SIM_TIME_STEPS             = 50;
LOG_NAME                            = ['k-Wave-Log-', getDateString];
HIGHEST_PRIME_FACTOR_WARNING        = 7;

% set additional literals that vary depending on the dimension of the
% simulation code
if exist('kgrid', 'var')
    switch kgrid.dim
        case 1
            PLOT_SCALE_DEF          = [-1.1, 1.1];
            PML_SIZE_DEF            = 20;
            PLOT_SCALE_WARNING      = 5;
        case 2
            PLOT_SCALE_DEF          = [-1, 1];
            PML_SIZE_DEF            = 20;
            PLOT_SCALE_WARNING      = 10;
        case 3
            PLOT_SCALE_DEF          = [-1, 1];
            PML_SIZE_DEF            = 10;
            PLOT_SCALE_WARNING      = 20;
    end
else
        PLOT_SCALE_DEF              = [-1, 1];
        PML_SIZE_DEF                = 20;
        PLOT_SCALE_WARNING          = 10;
end

% set the default CFL values used if kgrid.t_array is set to 'auto'
KSPACE_CFL                          = 0.3;
PSTD_CFL                            = 0.1;

% set additional literals unique to the elastic code
MULTI_AXIAL_PML_RATIO_DEF           = 0.1;

% =========================================================================
% HDF5 LIBRARY
% =========================================================================

% data type
DATA_TYPE_ATT_NAME                  = 'data_type';
MATRIX_DATA_TYPE_MATLAB             = 'single';
MATRIX_DATA_TYPE_C                  = 'float';
INTEGER_DATA_TYPE_MATLAB            = 'uint64';
INTEGER_DATA_TYPE_C                 = 'long';

% real / complex
DOMAIN_TYPE_ATT_NAME                = 'domain_type';
DOMAIN_TYPE_REAL                    = 'real';
DOMAIN_TYPE_COMPLEX                 = 'complex';

% file descriptors
FILE_MAJOR_VER_ATT_NAME             = 'major_version';
FILE_MINOR_VER_ATT_NAME             = 'minor_version';
FILE_DESCR_ATT_NAME                 = 'file_description';
FILE_CREATION_DATE_ATT_NAME         = 'creation_date';
CREATED_BY_ATT_NAME                 = 'created_by';

% file type
FILE_TYPE_ATT_NAME                  = 'file_type';
HDF_INPUT_FILE                      = 'input';
HDF_OUTPUT_FILE                     = 'output';
HDF_CHECKPOINT_FILE                 = 'checkpoint';

% file version information
HDF_FILE_MAJOR_VERSION              = '1';
HDF_FILE_MINOR_VERSION              = '1';

% compression level
HDF_COMPRESSION_LEVEL               = 0;