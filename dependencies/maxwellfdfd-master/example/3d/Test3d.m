inspect_only = false;
has_sol = false;  % true if solution files exist

% Solver Options
solveropts.method = 'inputfile';
filenamebase = 'tyler_3d';
solveropts.filenamebase = filenamebase;

if ~has_sol  % solution files do not exist
	% Input Files
	gray = [0.5 0.5 0.5];  % [r g b]
	[~, ~, obj_array, src_array] = maxwell_run(...
		'OSC', 1e-9, 1550, ...
		'DOM', {'Palik/SiO2', 'none'}, [-200 1700; -700 700; -600 600], 20, BC.p, 200, ...
		'OBJ', ...
			{'Palik/SiO2', 'none'}, Box([-200 1700; -50 50; -50 50], [20 2 2]), ...
			{'CRC/Ag', gray}, ...
				Box([-200 1700; -700 -25; -25 25], 20), Box([-200 1700; 25 700; -25 25], 20), ...
		'SRCJ', ModalSrc(Axis.x, 200), ...
		solveropts, inspect_only);

	if ~inspect_only
		save(filenamebase, 'obj_array', 'src_array');
	end
else  % solution files exist
	[E, H] = read_output(filenamebase);
	load(filenamebase);  % read obj_array and src_array
	visall(E{Axis.y}, obj_array, src_array);
end