
%% Solve the system.
inspect_only = false;
dl = 1;                         %grid resolution
hold on;
gray = [.5 .5 .5];

crange = (1:.1:1);
wrange = (2:1:2);               %spacing
hrange = (1:.1:1);               %width of dipole
srange = (150:50:150);
freqrange = (3e8/4000:1e4:3e8/200);
lamrange = 3e8./freqrange;
fprintf('Min Lambda   = %e\n', min(lamrange));
fprintf('Max Lambda   = %e\n', max(lamrange));
num_params = 5;                 % c, w, h, s vs lambda
num_coords = 3;
num_physical = 5;
total_steps = num_physical*length(freqrange)*length(srange)*length(hrange)*length(wrange)*length(crange);
steps = 0;

Estrength = zeros(num_physical, length(crange), length(wrange), length(hrange), length(srange), length(freqrange));
disp(size(Estrength));
Edir      = zeros(num_coords, length(crange), length(wrange), length(hrange), length(srange), length(freqrange));

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
                disp('working on 1');
                shape1 = [Box([-s-w/2, -w/2; -h/2, h/2; 0, 1]) Box([w/2, s+w/2; -h/2, h/2; 0, 1])];                 
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
                lami = 0;
                for freq = freqrange
                    lami = lami + 1;
                    lam = 3e8/freq;
                    steps = steps + 1;
                    disp(' ');
                    fprintf('SIZE   = %e\n', s);
                    fprintf('WIDTH  = %e\n', w);
                    fprintf('HEIGHT = %e\n', h);
                    fprintf('LAMBDA = %e\n', lam);
                    fprintf('C      = %e\n', c);
                    fprintf('PERCENTAGE DONE: %e\n', steps*100/total_steps);
                    disp(' ');
                    
                    %No branches
                    [E1, H1, ob_array1, src_array1, J1] = maxwell_run(...
                        'OSC', 1e-9, lam, ...
                        'DOM', {'vacuum', 'none', 1.0}, [-L-10 L+10; -L-10 L+10; 0 1], [dl,dl,1] , BC.p, [10 10 0], ...
                        'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 1], 1), ...
                        'SOBJ', ...  % scatter objects
                            {'perfect_conductor', gray, inf}, ...                               
                                shape5, ...
                            {'hole', gray, 1}, ...
                                hole, ...
                        'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 1], Axis.y, Axis.x), ...
                        inspect_only);
                    
                    steps = steps + 1;
                    disp(' ');
                    fprintf('SIZE   = %e\n', s);
                    fprintf('WIDTH  = %e\n', w);
                    fprintf('HEIGHT = %e\n', h);
                    fprintf('LAMBDA = %e\n', lam);
                    fprintf('C      = %e\n', c);
                    fprintf('PERCENTAGE DONE: %e\n', floor(steps*100/total_steps));
                    disp(' ');

                    [E2, H2, ob_array2, src_array2, J2] = maxwell_run(...
                        'OSC', 1e-9, lam, ...
                        'DOM', {'vacuum', 'none', 1.0}, [-L-10 L+10; -L-10 L+10; 0 1], [dl,dl,1] , BC.p, [10 10 0], ...
                        'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 1], 1), ...
                        'SOBJ', ...  % scatter objects
                            {'perfect_conductor', gray, inf}, ...
                                shape5, ...
                            {'hole', gray, 1}, ...
                                hole, ...
                                cover, ...
                        'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 1], Axis.y, Axis.x), ...
                        inspect_only);
                    
                    steps = steps + 1;
                    disp(' ');
                    fprintf('SIZE   = %e\n', s);
                    fprintf('WIDTH  = %e\n', w);
                    fprintf('HEIGHT = %e\n', h);
                    fprintf('LAMBDA = %e\n', lam);
                    fprintf('C      = %e\n', c);
                    fprintf('PERCENTAGE DONE: %e\n', floor(steps*100/total_steps));
                    disp(' ');
                    
                    %%
                                    %}
                    opts.withobjsrc = true;
                    %vis2d(E{Axis.x}, Axis.z, 0, ob_array, src_array);

                    [Xsize1,Ysize1,Zsize1] = size(E1{Axis.x}.array());
                    Xsize1 = floor(Xsize1/2);
                    Ysize1 = floor(Ysize1/2);
                    Zsize1 = floor(Zsize1/2);
                    
                    [Xsize2,Ysize2,Zsize2] = size(E2{Axis.x}.array());
                    Xsize2 = floor(Xsize2/2);
                    Ysize2 = floor(Ysize2/2);
                    Zsize2 = floor(Zsize2/2);                    

                    %}
                    Edir(1, 1, ci, wi, hi, si, lami) = abs(E1{Axis.x}.array(Xsize1,Ysize1,Zsize1))^2;
                    Edir(1, 2, ci, wi, hi, si, lami) = abs(E1{Axis.y}.array(Xsize1,Ysize1,Zsize1))^2;
                    Edir(1, 3, ci, wi, hi, si, lami) = abs(E1{Axis.z}.array(Xsize1,Ysize1,Zsize1))^2;
                    
                    Edir(2, 1, ci, wi, hi, si, lami) = abs(E2{Axis.x}.array(Xsize2,Ysize2,Zsize2))^2;
                    Edir(2, 2, ci, wi, hi, si, lami) = abs(E2{Axis.y}.array(Xsize2,Ysize2,Zsize2))^2;
                    Edir(2, 3, ci, wi, hi, si, lami) = abs(E2{Axis.z}.array(Xsize2,Ysize2,Zsize2))^2;

                    disp(wi);
                    %}
                    Estrength(1, ci, wi, hi, si, lami) = ...
                                      abs(E1{Axis.x}.array(Xsize1,Ysize1,Zsize1))^2 ...
                                    + abs(E1{Axis.y}.array(Xsize1,Ysize1,Zsize1))^2 ...
                                    + abs(E1{Axis.z}.array(Xsize1,Ysize1,Zsize1))^2;
                    disp(Estrength(1, ci, wi, hi, si, lami));
                    disp(ci);
                    disp(wi);
                    disp(hi); 
                    disp(si);
                    disp(lami);
                                
                    Estrength(2, ci, wi, hi, si, lami) = ...
                                      abs(E2{Axis.x}.array(Xsize2,Ysize2,Zsize2))^2 ...
                                    + abs(E2{Axis.y}.array(Xsize2,Ysize2,Zsize2))^2 ...
                                    + abs(E2{Axis.z}.array(Xsize2,Ysize2,Zsize2))^2;

                end
            end
        end
    end
end
%%
%lines = ['k--.'; 'g--.'; 'b--.'; 'r--.'; 'k-.'; 'g-.'; 'b-.'; 'r-.'];
figure(10);
hold on;
title('Branches');
xlabel('wavelength (nm)');
ylabel('|E(0,0,0)|^2');
l = [];

if (num_physical > 1)
    for i = (1:num_physical)
        plot(3./lamrange, squeeze(Estrength(i,1,1,1,1,:)));
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
