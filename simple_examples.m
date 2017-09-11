clear all; clc;
% This is intended to help a new user get started with the simulation
% class.  
% To get started working on this, you must be inside of the main directory.
% For a more detailed view of the methods and properties of this
% class, please type 'doc simulation' into your matlab command line.

%% This is the simplest simulation you could run.

s = simulation;             % Begin the simulation (SOI by default)
% s = simulation('SiN')     % Begin the simulation with Si3N4/SiO2 materials
% s = simulation('SOI')     % Begin the simulation with SOI materials

s.setup;                    % solve waveguides, geometry, create input pulse, solve for DLA response
s.run_simulation;           % do input coupling, propagate pulse, compute gradient   

display(['final gradient of ', num2str(s.G/1e6), ' MV/m.']);

%% Here is the above simulation broken down into its subcomponents 
% (equivalent in practice but if changing only one of the parameters at a
% time, may make more sense to simly use one of these submethods, assuming
% the others have been run already

s = simulation;

% The methods run by s.setup
s.solve_waveguide_parameters;               % takes wg geometry and material definitions, computes group index, dispersion parameters for NLSE and geometry solver {s.n_eff, s.n_g, s.A_eff, s.beta2}
s.solve_geometry;                           % takes group index and beta and figures out a structure to velocity match pulses to e- beam {s.Ls, s.Rs, s.thetas} 

% the methods run by s.run_simulation
s.compute_P_max;                            % compute minimum fields before constraints are reached
s.compute_gradient;                         % computes the gradient at the end of the structure using lorenzian fit of graident(f) {s.G, s.G_scan, s.E, s.damaged}

% runs all of the plotting methods
close all;
s.plot_waveguide;                           % plots waveguide mode, dispersion curves.
s.plot_damage;                            % plots geometry

%% When changing any releavant parameter away from default, you must re-run some (or all) of the setup methods.  For example

s.tau = 300e-15;              % new pulse duration of 300 fs.  (default = 250 fs)
s.wg_height = 250e-9;         % new waveguide core height of 250 nm.  (default = 200 nm for SOI)
s.solve_waveguide_parameters; % resolve for the waveguide parameters.
s.L = 96e-6;                  % new interaction length of 96 um
s.solve_geometry;             % new group index, so need to re-solve structure.


s.run_simulation;             % after setup, we can run_simulation again
display(['final gradient of ', num2str(s.G/1e6), ' MV/m.']);


