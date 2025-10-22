function setup_oceanmesh2d()
% add paths for OceanMesh2D

BasePath = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(BasePath,'datasets')))
addpath(genpath(fullfile(BasePath,'utilities')))
