#!/bin/sh
# Initial job file, fill out for for ultrasound simulations

#PBS -lwalltime=01:01:10

#PBS -lselect=1:ncpus=16:mem=96gb:ngpus=4:gpu_type=RTX6000

module load matlab

matlab -nodisplay -nosplash -nodesktop -r "run('/rds/general/user/sk8717/home/UltrasoundProject/simulations/main_files/Main_simulation_3D.m')"
