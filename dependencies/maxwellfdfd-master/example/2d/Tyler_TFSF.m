clear all; close all; clear classes; clc;

%% Set flags.
inspect_only = false;

%% Solve the system.
w = 25;
dL = 10;
dl = 2;
h = 40;
wvlen = 200;
clear solveropts;
hold off;
hold on;
ps = [];
i = [];
ones = [];
legends = [];
bs = [];
srange = (100:10:100);
for s = srange
    lams = [];
    Estrength = [];
    Estrength2 = [];
    p2 = [];
    b = (s + w/2)*2;
    bs = [bs, b];
    for lam = (400:5:1300)
        fprintf('SIZE   = %e\n', s*2);
        fprintf('LAMBDA = %e\n', lam);
        [E, H, ob_array, src_array, J] = maxwell_run(...
            'OSC', 1e-9, lam, ...
            'DOM', {'vacuum', 'none', 1.0}, [-610 610; -610 610; 0 dl], [dL dL dl], BC.p, [10*dL 10*dL 0], ...
            'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 dl], dl), ...
            'SOBJ', ...  % scatter objects
                {'Johnson/Au', 'y'}, ...
                    PolygonalCylinder(Axis.z, dl, dl/2, [w/2 0; w/2+s*sqrt(3)/2 -s/2; w/2+s*sqrt(3)/2 s/2], dl), ...
                    PolygonalCylinder(Axis.z, dl, dl/2, [-w/2 0; -w/2-s*sqrt(3)/2 s/2; -w/2-s*sqrt(3)/2 -s/2], dl), ...
            'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 dl], Axis.y, Axis.x), ...
            inspect_only);
        [E2, H2, ob_array2, src_array2, J2] = maxwell_run(...
            'OSC', 1e-9, lam, ...
            'DOM', {'vacuum', 'none', 1.0}, [-610 610; -610 610; 0 1], [2,2,1] , BC.p, [10 10 0], ...
            'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 1], 1), ...
            'SOBJ', ...  % scatter objects
                {'Johnson/Au', 'y'}, ...
                    Box([-s-w/2, -w/2; -h/2, h/2; 0, 1]), Box([w/2, s+w/2; -h/2, h/2; 0, 1]), ...  % metal slit
            'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 1], Axis.y, Axis.x), ...
            inspect_only);
  %      [E2, H2, obj_array2, src_array2, J2] = maxwell_run(...
  %          'OSC', 1e-9, lam, ...
  %          'DOM', {'vacuum', 'none', 1.0}, [-610 610; -610 610; 0 dl], [dL dL dl], BC.p, [10*dL 10*dL 0], ...
  %          'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 dl], dl), ...
  %          'SOBJ', ...  % scatter objects
  %              {'vacuum', 'none', 1.0}, ...
  %                  PolygonalCylinder(Axis.z, dl, dl/2, [w/2 0; w/2+s*sqrt(3)/2 -s/2; w/2+s*sqrt(3)/2 s/2], dl), ...
  %                  PolygonalCylinder(Axis.z, dl, dl/2, [-w/2 0; -w/2-s*sqrt(3)/2 s/2; -w/2-s*sqrt(3)/2 -s/2], dl), ...
  %          'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 dl], Axis.y, Axis.x), ...
  %          inspect_only);
        %% Visualize the solution.
        %figure
        %clear opts
        %opts.withobjsrc = true;
        %opts.withabs = true;
        % opts.withinterp = false;
        % opts.withgrid = true;
        % opts.cscale = 1e-1;
        % opts.cmax = 1.4;
        %z_location = 0;
        % vis2d(E{Axis.x}, Axis.z, z_location, obj_array, src_array, opts)
        % vis2d(H{Axis.z}, Axis.z, z_location, obj_array, src_array, opts)

        % %% Calculate the power emanating from the source.
        power1 = powerflux_box(E,H,[-10 10; -10 10; 0 1]);
  %      power2 = powerflux_box(E2,H2,[-10 10; -10 10; 0 1]);
        %Estrength = [Estrength, abs(E{Axis.x}.array(141,103,1))^2 ...
        %        + abs(E{Axis.y}.array(141,103,1))^2 ...
        %        + abs(E{Axis.z}.array(141,103,1))^2  ];
        Size = size(E{Axis.x}.array());
        Estrength = [Estrength, abs(E{Axis.x}.array(floor(Size(1)/2),floor(Size(2)/2),1))^2 ];
        Estrength2 = [Estrength2, abs(E2{Axis.x}.array(floor(Size(1)/2),floor(Size(2)/2),1))^2 ];        
        lams = [lams,lam];
      %  p1 = [p1,power1];
  %      p2 = [p2,power2];
         %fprintf('power = %e\n', power);
    end
%    legends = [legends, num2str((s+w/2)*2)];
%    ones = [ones, 1];
    ps = [ps, [Estrength]];
    [M,I] = max(Estrength);
    [M2,I2] = max(Estrength2);
    i = [i,I];
    hold on
    plot(lams, Estrength);
    plot(lams, Estrength1);
end
%plot(srange*2, lams(i));
[slope, intercept] = polyfit(srange*2, lams(i), 1);
%stem(legends,ones);
