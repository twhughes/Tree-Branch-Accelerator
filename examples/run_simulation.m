%% Simulation example script - Tyler Hughes
% use 'doc simulation' for a cleaner view of the methods and properties if the simulation class 

% ---CLEAR STUFF AND LOAD WHOLE DIRECTORY---%
clear all; clc;
addpath(genpath('.'))

% ----INITIALIZE----%

% initialize a simulation object
s = simulation;
% if s.verbose, then it will print out status and calculations (good for
% debugging, annoying for multiple loops).
s.verbose = false; 
s.verbose = true;

% ----WAVEGUIDE SETUP----%

% specify material stack and corresponding dimensions ('SOI' or 'SiN')
s.material_stack = 'SOI';
% material_stack = 'SiN';
% load in the wavegudie dimensions based on which material stack you want
if strcmp(s.material_stack,'SOI')
    s.wg_width = 1e-6;
    s.wg_height = 200e-9;
    s.clad_height = 600e-9;
else strcmp(s.material_stack,'SiN')
    s.wg_width = 4e-6;
    s.wg_height = 100e-9;
    s.clad_height = 5e-6;
end

% if any of the other parameters (n_eff, n_g, beta2, A_eff) are known,
% one may pre-specify them to save time.
% If all are not defined, we will solve with the waveguide mode solver.
s.solve_waveguide_parameters();

%% GEOMETRY
L = 1e-3;       % total end interaction length (m)
beta = 1;       % electron speed/speed of light
M = 2;          % number of DLA periods powered by 1 waveguide
L0 = 1e-3;      % distance from input coupler to first split (m)

% solve for the resulting lengths, radii, and angles
s.solve_geometry();
%s.draw_geometry();             % uncomment if you wish to draw the geometry

%% SPLITTING

% directly set constant split power efficiency
s.split_efficiency = 0.9;

%% BEND RADIUS

% Need to DO! assumed bends to not give loss for now

%% Damage

% using waveguide materials, load damage threshold as fn. of pulse duration
s.load_damage();                        % 's.Ed(tau)' gives the damage threshold (V/M) as a fn. of tau (s)

%% PROPAGATE PULSE
% input pulse parameters
lambda = 2e-6;
E0 = 1e10;
tau = 250e-15;
tbp = 0.44;
chirp = 0;
rep_rate = 1e6;

% construct A0_t (initial pulse)
s.construct_initial_pulse();
s.input_couple();
% propagate pulse through structure.  The end pulse is stored in s.A_t 
s.propagate_pulse();

%% DLA section

% find the pillar radius to maximize gradient at central frequency
% then do frequenct scan to find Q factor and enhancement
s.NR = 100;
s.NF = 100;
s.R_pillars = 1.597e-6;
s.Q = 173.866;
s.T0 = 0.63442;
s.compute_DLA_parameters();
s.compute_gradient();

display(['gradient    = ', num2str(s.G/1e6), ' MV/m']);
display(['energy gain = ', num2str(s.E_gain/1e3), ' keV']);










