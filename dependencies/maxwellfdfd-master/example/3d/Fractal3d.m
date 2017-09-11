inspect_only = false;
has_sol = true;  % true if solution files exist

% Solver Options
solveropts.method = 'inputfile';

h = 2;
w = 2;
s = 50;
b = 3*s;
c = 1;
lamrange = (500:500);
z = 6;%max(lamrange);

%generate fractal geometries
global A;
A = [];
shape = [Box([-s-w/2, -w/2; -h/2, h/2; 0, 1]) Box([w/2, s+w/2; -h/2, h/2; 0, 1])];                 
A = [];
disp('working on 2');                
HTree2(1, s, h, c);
shape2 = A;
A = [];
disp('working on 3');
HTree2(2, s, h, c);
shape3 = A;      
A = [];
disp('working on 4');
HTree2(3, s, h, c);
shape4 = A;  
A = [];
disp('working on 5');
HTree2(4, s, h, c);
shape5 = A;

gray = [0.5 0.5 0.5];  % [r g b]
%ModalSrc(Axis.z, -b*2/3), ...
for lam = lamrange
    
    filenamebase = strcat('fractal_3d_', int2str(lam));
    solveropts.filenamebase = filenamebase;
    if ~has_sol  % solution files do not exist
        % Input Files
        [~, ~, obj_array, src_array] = maxwell_run(...
            'OSC', 1e-9, lam, ...
            'DOM', {'Air', 'none', 1+0j}, [-b b; -b b; -z z], 1, BC.p, 2, ...
            'SOBJ', ...
                {'Air', 'none', 1+0j}, Box([-b*2/3 b*2/3; -b*2/3 b*2/3; -z*2/3 z*2/3], [1 1 1]), ...
                {'Perfect Conductor', gray, inf}, ...
                    shape3, ...
            'SRCJ', TFSFPlaneSrc([-b b; -b b; -z z], Axis.z, Axis.x), ...
            solveropts, inspect_only);

        if ~inspect_only
            save(filenamebase, 'obj_array', 'src_array');
        end
    else  % solution files exist
        [E, H] = read_output(filenamebase);
        load(filenamebase);  % read obj_array and src_array
        visall(E{Axis.y}, obj_array, src_array);
    end
end