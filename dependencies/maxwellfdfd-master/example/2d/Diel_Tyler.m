%% Set flags.
inspect_only = false;
b = 500; %inner box L/2
c = 10;
t = 150; %thickness
p_up = [];
p_down = [];
p_up2 = [];
p_down2 = [];
Pp_up = [];
Pp_down = [];
Pp_up2 = [];
Pp_down2 = [];
freqs = (.05:.005:1.0);
lams = fliplr(1./freqs);
%% Solve the system.
for lam = lams
    gray = [0.5 0.5 0.5];  % [r g b]
    [E, H, obj_array, src_array, J] = maxwell_run(...
        'OSC', 1e-9, lam, ...
        'DOM', {'vacuum', 'none', 1}, [0 1; -2*b 2*b; 0 1], [1,5,1], BC.p, [0 b/5 0], ...
        'OBJ', {'vacuum', 'none', 1.0}, Box([0 1; -2*b 2*b; 0 1], 1), ...
        'SOBJ', ...  % scatter objects
            {'vacuum' 'none', 4}, Box([0 1; -t/2 t/2; 0 1]), ...
           'SRCJ', TFSFPlaneSrc([0 1; -b b; 0 1], Axis.y, Axis.z), ...
        inspect_only);

    
   [E2, H2, obj_array2, src_array2, J2] = maxwell_run(...
        'OSC', 1e-9, lam, ...
        'DOM', {'vacuum', 'none', 1}, [0 1; -2*b 2*b; 0 1], [1,5,1], BC.p, [0 b/5 0], ...
        'OBJ', ...  % scatter objects
            {'vacuum' 'none', 4}, Box([0 1; -t/2 t/2; 0 1]), ...
           'SRCJ', PlaneSrc(Axis.y, -b, Axis.z), ...
        inspect_only);
    
    
    %% Visualize the solution.
    if ~inspect_only
    % 	opts.withgrid = true;
        opts.withobjsrc = true;
        opts.withabs = true;
        opts.withpml = false;
        opts.phase = pi/2;
      %  vis2d(E{Axis.x}, Axis.z, 0.5, obj_array, src_array, opts)

        p_up = [p_up, powerflux_patch(E, H, Axis.y, b/2)];
        p_down = [p_down, powerflux_patch(E, H, Axis.y, -b/2)];
        p_up2 = [p_up2, powerflux_patch(E, H, Axis.y, 3*b/2)];
        p_down2 = [p_down2, powerflux_patch(E, H, Axis.y, -3*b/2)];
        
        Pp_up = [Pp_up, powerflux_patch(E2, H2, Axis.y, b/2)];
        Pp_down = [Pp_down, powerflux_patch(E2, H2, Axis.y, -b/2)];
        Pp_up2 = [Pp_up2, powerflux_patch(E2, H2, Axis.y, 3*b/2)];
        Pp_down2 = [Pp_down2, powerflux_patch(E2, H2, Axis.y, -3*b/2)];
    end
end
%%
c0 = 299000000;
hold on
plot(c0./lams, p_up, 'ro-');
plot(c0./lams, p_down, 'bo-');
plot(c0./lams, p_up2, 'go-');
plot(c0./lams, p_down2, 'ko-');
plot(c0./lams, Pp_up, 'r-');
plot(c0./lams, Pp_down, 'b-');
plot(c0./lams, Pp_up2, 'g-');
plot(c0./lams, Pp_down2, 'k-');