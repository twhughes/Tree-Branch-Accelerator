clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

%% Create shapes.
a = 100;  % lattice constant
t = 1;  % slab thickness
r = 0.29*a;  % hole radius

ad = 25;  % divider for a
td = 10;  % divider for t
dd = 10;  % divider for d = 2*r

mx = 20.5;
my = 7.5;
slab_yn = Box([-mx*a mx*a; -my*a -0.5*a; 0 t], [a/ad, a/ad, t]);
slab_yp = Box([-mx*a mx*a; 0.5*a my*a; 0 t], [a/ad, a/ad, t]);

rod = CircularCylinder(Axis.z, t, [0 0 t/2], r, [2*r/dd, 2*r/dd, t]);
lams = (50:5:600);
pup = [];
pdn = [];
%% Solve the system.
gray = [0.5 0.5 0.5];  % [r g b]
src_loc = -5*a;
for lam = lams
    [E, H, obj_array, src_array, J] = maxwell_run(...
        'OSC', 1e-9, lam, ...
        'DOM', {'vacuum', 'none', 1}, [-10 10; -1000 1000; 0 t], [10 10 t], BC.p, [2 2 0], ...
        'OBJ', ...
            {'diel', gray, 2}, Box([-10 10; -a/2 a/2; 0 t]), ...
        'SRCJ', PlaneSrc(Axis.y, -500 , Axis.z), ...
            inspect_only);

    %% Visualize the solution.
    if ~inspect_only
        clear opts
    % 	opts.withgrid = true;
        opts.withobjsrc = true;
        opts.withabs = true;
        opts.withpml = false;
        opts.phase = pi/2;
    %    vis2d(E{Axis.z}, Axis.z, 0.5, obj_array, src_array, opts)
        %%
    % 	figure(2);
    % 	vis2d(H{Axis.y}, Axis.z, 0.5, obj_array, src_array, opts)

        %%
        flux_loc = 3*a;
        power_right = powerflux_patch(E, H, Axis.y, 500);
        power_left = -powerflux_patch(E, H, Axis.y, -500);
        fprintf('power:\n');
        fprintf('right = %s\n', num2str(power_right));
        fprintf('left = %s\n', num2str(power_left));
        fprintf('error = %s%%\n',num2str((power_left-power_right)/power_right*100));
        pup = [pup, power_right];
        pdn = [pdn, power_left];
        %%
    %	Sx = poynting(Axis.x, E{Axis.y}, E{Axis.z}, H{Axis.y}, H{Axis.z}, Axis.y, 0);
    %	[array, l] = Sx.data_original;
        %figure(3);
        %plot(l{2}, abs(array))
        mx*a - 10*a;
    end
end
figure(3);
hold on;
plot(lams,pdn);
plot(lams,pup);