
s = simulation_simple;

s.material_stack = 'SOI';
s.wg_width = 1e-6;
s.wg_height = 200e-9;
s.clad_height = 600e-9;

s.M = 4;
s.L = s.M*s.lambda*s.beta*2^5;
s.L
s.split_efficiency = 1;%0.95;
s.bend_efficiency = 1;%0.95;
s.coupler_efficiency = 1;
s.phi_SPM_limit = 2*pi;

s.num_taus_sample = 500;
s.NT = 2^15;
s.solve_waveguide_parameters;
%close all;
%s.plot_waveguide;
%s.solve_geometry;

Q_tau;
%s.tau = 10^1.5*1e-12;
%s.Q = 10^1;
%s.compute_P_max;
%s.compute_gradient;

%%
gradient_MV_m = max(max(abs(Gs)))/1e6
Energy_keV = max(max(abs(E_gains)))/1e3
num_needed = 1e6/max(max(abs(E_gains)))
length_um = s.L*1e6