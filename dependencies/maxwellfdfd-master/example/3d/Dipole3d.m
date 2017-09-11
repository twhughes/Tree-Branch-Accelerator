clear all; close all; clear classes; clc;

%% Set flags.
is_new = true;
inspect_only = false;
filenamebase = 'slot_3d_basic';

if is_new
	%% Create input files.
	gray = [0.5 0.5 0.5];  % [r g b]
	solveropts.method = 'inputfile';
	solveropts.filenamebase = filenamebase;
	solveropts.showstruct = false;
    b = 600; %domain size
	% (along x)
	modeopts.clue = 'guess';
	modeopts.neff = 2;
	[E, H, obj_array, src_array, J] = maxwell_run(...
		'OSC', 1e-9, 1550, ...
		'DOM', {'air', 'none', 0}, [-b b; -b b; -b b], 5, BC.p, b/5, ...
		'OBJ', ...
			{'air', 'none', 0}, Box([-b/2 b/2; -b/2 b/2; -b/2 b/2], [20 20 20]), ...
			{'Perfect Conductor', gray, inf}, ...
				Box([-200 -10; -10 10; -10 10], 5), Box([10 200; -10 10; -10 10], 5), ...
		'SRCJ', ModalSrc(Axis.z, 300, modeopts), ...
		solveropts, inspect_only);
	
	save(filenamebase, 'obj_array', 'src_array');
else
	[E, H] = read_output(filenamebase);
	load(filenamebase);
end

%% Visualize the solution.
if ~is_new
	figure;
	clear opts;
	opts.withabs = true;
	opts.isopaque = false;
% 	opts.withgrid = true;
	visall(E{Axis.x}, obj_array, src_array, opts);
end
