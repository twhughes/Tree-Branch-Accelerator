
%% Solve the system.
inspect_only = false;
dl = 1;                         %grid resolution
hold on;
gray = [.5 .5 .5];

crange = (1:.1:1);
wrange = (2:1:2);               %spacing
hrange = (2:2:2);               %width of dipole
srange = (50:50:50);
lamrange = (200:100:400);
num_params = 5;                 % c, w, h, s vs lambda
num_coords = 3;
num_physical = 2;
total_steps = num_physical*length(lamrange)*length(srange)*length(hrange)*length(wrange)*length(crange);
steps = 0;

Estrength = zeros(num_physical, length(crange), length(wrange), length(hrange), length(srange), length(lamrange));
disp(size(Estrength));
Edir      = zeros(num_coords, length(crange), length(wrange), length(hrange), length(srange), length(lamrange));
shape = [];

si   = 0;
for s = srange
    si = si + 1;
    wi   = 0;
    for w = wrange
        wi = wi + 1;
        b = (s + w/2)*2;
        L = b*3/2;
        hi   = 0;
        for h = hrange
            hi = hi + 1;
            hole = Box([-w/2 w/2; -h h; 0 1]);
            cover = Box([-b+h b-h; -b+h -h/2; 0 1]);            
            ci   = 0;
            for c = crange
                ci = ci + 1;
               	global A;
                A = [];
                shape1 = [Box([-s-w/2, -w/2; -h/2, h/2; 0, 1]), Box([w/2, s+w/2; -h/2, h/2; 0, 1])];              
                for p = (2:num_physical+1)                 
                    A = [];
                    HTree2(p, s, h, c);
                    shape = [shape, A];
                end
                lami = 0;
                for lam = lamrange
                    lami = lami + 1;
                    
                    steps = steps + 1;
                    disp(' ');
                    fprintf('SIZE   = %e\n', s);
                    fprintf('WIDTH  = %e\n', w);
                    fprintf('HEIGHT = %e\n', h);
                    fprintf('LAMBDA = %e\n', lam);
                    fprintf('C      = %e\n', c);
                    fprintf('PERCENTAGE DONE: %e\n', steps*100/total_steps);
                    disp(' ');
                       %initialize field arrays
                    E         = zeros(num_physical);
                    H         = zeros(num_physical);
                    ob_array  = zeros(num_physical);
                    src_array = zeros(num_physical);
                    J         = zeros(num_physical);
                    
                    for p = (1:num_physical)
                        %No branches
                        if p == 1
                            myshape = shape1;
                        else
                            myshape = shape(p);
                        end
                         
                        [Etemp, Htemp, ob_arraytemp, src_arraytemp, Jtemp] = maxwell_run(...
                            'OSC', 1e-9, lam, ...
                            'DOM', {'vacuum', 'none', 1.0}, [-L-10 L+10; -L-10 L+10; 0 1], [dl,dl,1] , BC.p, [10 10 0], ...
                            'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 1], 1), ...
                            'SOBJ', ...  % scatter objects
                                {'perfect_conductor', gray, inf}, ...
                                    myshape, ...
                                {'hole', 'none', 1}, ...
                                    hole, ...
                                    cover, ...
                            'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 1], Axis.y, Axis.x), ...
                            inspect_only);

                        steps = steps + 1;
                        disp(' ');
                        fprintf('PHYS   = %e\n', p);                        
                        fprintf('SIZE   = %e\n', s);
                        fprintf('WIDTH  = %e\n', w);
                        fprintf('HEIGHT = %e\n', h);
                        fprintf('LAMBDA = %e\n', lam);
                        fprintf('C      = %e\n', c);
                        fprintf('PERCENTAGE DONE: %e\n', floor(steps*100/total_steps));
                        disp(' ');

                        %%
                        
                        opts.withobjsrc = true;

                        [Xsize,Ysize,Zsize] = size(Etemp{Axis.x}.array());
                        Xsize = floor(Xsize/2);
                        Ysize = floor(Ysize/2);
                        Zsize = floor(Zsize/2);
 
                        for q = (1:num_coords)
                           Edir(p, q, ci, wi, hi, si, lami) = abs(Etemp{Axis.x}.array(Xsize,Ysize,Zsize))^2;
                        end
                        %}
                        vis2d(Etemp{Axis.x}, Axis.z, 0, ob_arraytemp, src_arraytemp);
                        Estrength(p, ci, wi, hi, si, lami) = ...
                                          abs(Etemp{Axis.x}.array(Xsize,Ysize,Zsize))^2 ...
                                        + abs(Etemp{Axis.y}.array(Xsize,Ysize,Zsize))^2 ...
                                        + abs(Etemp{Axis.z}.array(Xsize,Ysize,Zsize))^2;

                    end
                end
            end
        end
    end
end
%%
%lines = ['k--.'; 'g--.'; 'b--.'; 'r--.'; 'k-.'; 'g-.'; 'b-.'; 'r-.'];
figure(1);
hold on;
title('Branches');
xlabel('wavelength (nm)');
ylabel('|E(0,0,0)|^2');
l = [];

if (num_physical > 1)
    for i = (1:num_physical)
        plot(lamrange, squeeze(Estrength(i,1,1,1,1,:)));
        l = [legend, num2str(2^i)];
    end
end

if (ci > 1)
    figure(2);
    hold on;
    title('C');
    xlabel('wavelength (nm)');
    ylabel('|E(0,0,0)|^2');
    l = num2str(crange);
    for i = (1:ci)
        plot(lamrange, squeeze(Estrength(1,i,1,1,1,:)));
    end
end

if (wi > 1)
    figure(3);
    hold on;
    title('Width');
    xlabel('wavelength (nm)');
    ylabel('|E(0,0,0)|^2');
    l = num2str(wrange);
    for i = (1:wi)
        plot(lamrange, squeeze(Estrength(5,1,i,1,1,:)));
    end
end
   
if (hi > 1)
    figure(4);
    hold on;
    title('Height');
    xlabel('wavelength (nm)');
    ylabel('|E(0,0,0)|^2');
    l = num2str(hrange);
    for i = (1:hi)
        plot(lamrange, squeeze(Estrength(1,1,1,i,1,:)));
    end
end

if (si > 1)
    figure(5);
    hold on;
    title('Length');
    xlabel('wavelength (nm)');
    ylabel('|E(0,0,0)|^2');
    l = num2str(srange);
    for i = (1:si)
        plot(lamrange, squeeze(Estrength(1,1,1,1,i,:)));
    end
end
