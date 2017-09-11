
%% Solve the system.
inspect_only = false;
dl = 1;                         %grid resolution
hold on;
gray = [.5 .5 .5];

crange = (1:.1:1);
wrange = (1:1:3);               %spacing
hrange = (1:.1:1);               %width of dipole
srange = (100:10:100);
lamrange = (200:50:1800);
num_params = 5;                 % c, w, h, s vs lambda
num_coords = 3;
num_physical = 5;
total_steps = num_physical*length(lamrange)*length(srange)*length(hrange)*length(wrange)*length(crange);
steps = 0;

Estrength = zeros(num_physical, length(crange), length(wrange), length(hrange), length(srange), length(lamrange));
disp(size(Estrength));
Edir      = zeros(num_coords, length(crange), length(wrange), length(hrange), length(srange), length(lamrange));
global A;
paramsList = ndgrid(crange, wrange, hrange, srange, lamrange);
disp(paramsList);
parfor s = srange
    for w = wrange
        b = (s + w/2)*2;
        L = b*3/2;
        for h = hrange
            hole = Box([-w/2 w/2; -h h; 0 1]);
            for c = crange
                A = [];
                shape1 = [Box([-s-w/2, -w/2; -h/2, h/2; 0, 1]) Box([w/2, s+w/2; -h/2, h/2; 0, 1])];                 
                A = [];
                HTree2(1, s, h, c);
                shape2 = A;
                A = [];
                HTree2(2, s, h, c);
                shape3 = A;      
                A = [];
                HTree2(3, s, h, c);
                shape4 = A;  
                A = [];
                HTree2(4, s, h, c);
                shape5 = A;
           %     lami = 0;
                for lam = lamrange
             %       lami = lami + 1;
                    
              %      steps = steps + 1;
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
                                shape1, ...
                            {'hole', 'none', 1}, ...
                                hole, ...
                        'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 1], Axis.y, Axis.x), ...
                        inspect_only);
                    
                    disp( abs(E1{Axis.x}.array(Xsize1,Ysize1,Zsize1))^2 ...
                                    + abs(E1{Axis.y}.array(Xsize1,Ysize1,Zsize1))^2 ...
                                    + abs(E1{Axis.z}.array(Xsize1,Ysize1,Zsize1))^2)
                                
                 %   steps = steps + 1;
                    disp(' ');
                    fprintf('SIZE   = %e\n', s);
                    fprintf('WIDTH  = %e\n', w);
                    fprintf('HEIGHT = %e\n', h);
                    fprintf('LAMBDA = %e\n', lam);
                    fprintf('C      = %e\n', c);
                    fprintf('PERCENTAGE DONE: %e\n', floor(steps*100/total_steps));
                    disp(' ');
                    
                    %2^1 Branches
                    [E2, H2, ob_array2, src_array2, J2] = maxwell_run(...
                        'OSC', 1e-9, lam, ...
                        'DOM', {'vacuum', 'none', 1.0}, [-L-10 L+10; -L-10 L+10; 0 1], [dl,dl,1] , BC.p, [10 10 0], ...
                        'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 1], 1), ...
                        'SOBJ', ...  % scatter objects
                            {'perfect_conductor', gray, inf}, ...
                                shape2, ...
                            {'hole', 'none', 1}, ...
                                hole, ...
                            {'perfect_conductor', gray, inf}, ...
                        'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 1], Axis.y, Axis.x), ...
                        inspect_only);
                    
           %         steps = steps + 1;
                    disp(' ');
          
                    fprintf('SIZE   = %e\n', s);
                    fprintf('WIDTH  = %e\n', w);
                    fprintf('HEIGHT = %e\n', h);
                    fprintf('LAMBDA = %e\n', lam);
                    fprintf('C      = %e\n', c);
                    fprintf('PERCENTAGE DONE: %e\n', floor(steps*100/total_steps));
                    disp(' ');
                                      
                    %2^2 Branches
                    [E3, H3, ob_array3, src_array3, J3] = maxwell_run(...
                        'OSC', 1e-9, lam, ...
                        'DOM', {'vacuum', 'none', 1.0}, [-L-10 L+10; -L-10 L+10; 0 1], [dl,dl,1] , BC.p, [10 10 0], ...
                        'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 1], 1), ...
                        'SOBJ', ...  % scatter objects
                            {'perfect_conductor', gray, inf}, ...
                                shape3, ...
                            {'hole', 'none', 1}, ...
                                hole, ...
                        'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 1], Axis.y, Axis.x), ...
                        inspect_only);        

           %         steps = steps + 1;
                    disp(' ');
                    fprintf('SIZE   = %e\n', s);
                    fprintf('WIDTH  = %e\n', w);
                    fprintf('HEIGHT = %e\n', h);
                    fprintf('LAMBDA = %e\n', lam);
                    fprintf('C      = %e\n', c);
                    fprintf('PERCENTAGE DONE: %e\n', floor(steps*100/total_steps));
                    disp(' ');
                                     
                    %2^3 Branches
                    [E4, H4, ob_array4, src_array4, J4] = maxwell_run(...
                        'OSC', 1e-9, lam, ...
                        'DOM', {'vacuum', 'none', 1.0}, [-L-10 L+10; -L-10 L+10; 0 1], [dl,dl,1] , BC.p, [10 10 0], ...
                        'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 1], 1), ...
                        'SOBJ', ...  % scatter objects
                            {'perfect_conductor', gray, inf}, ...
                                shape4, ...
                            {'hole', 'none', 1}, ...
                                hole, ...
                        'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 1], Axis.y, Axis.x), ...
                        inspect_only);  
                    
                    shape5 = [Box([-s-w/2, -w/2; -h/2, h/2; 0, 1]) Box([w/2, s+w/2; -h/2, h/2; 0, 1])];
                    [E5, H5, ob_array5, src_array5, J5] = maxwell_run(...
                        'OSC', 1e-9, lam, ...
                        'DOM', {'vacuum', 'none', 1.0}, [-L-10 L+10; -L-10 L+10; 0 1], [dl,dl,1] , BC.p, [10 10 0], ...
                        'OBJ', {'vacuum', 'none', 1.0}, Box([-w/2 w/2; -w/2 w/2; 0 1], 1), ...
                        'SOBJ', ...  % scatter objects
                            {'none', gray, 1}, ...
                                shape1, ...
                            {'hole', 'none', 1}, ...
                                hole, ...
                        'SRCJ', TFSFPlaneSrc([-b b; -b b; 0 1], Axis.y, Axis.x), ...
                        inspect_only);   
                    %%
                                    %}
               %     opts.withobjsrc = true;
                    %vis2d(E{Axis.x}, Axis.z, 0, ob_array, src_array);
%{
                    [Xsize1,Ysize1,Zsize1] = size(E1{Axis.x}.array());
                    Xsize1 = floor(Xsize1/2);
                    Ysize1 = floor(Ysize1/2);
                    Zsize1 = floor(Zsize1/2);
                    
                    [Xsize2,Ysize2,Zsize2] = size(E2{Axis.x}.array());
                    Xsize2 = floor(Xsize2/2);
                    Ysize2 = floor(Ysize2/2);
                    Zsize2 = floor(Zsize2/2);                    
                    [Xsize3,Ysize3,Zsize3] = size(E3{Axis.x}.array());
                    Xsize3 = floor(Xsize3/2);
                    Ysize3 = floor(Ysize3/2);
                    Zsize3 = floor(Zsize3/2);
                    [Xsize4,Ysize4,Zsize4] = size(E4{Axis.x}.array());
                    Xsize4 = floor(Xsize4/2);
                    Ysize4 = floor(Ysize4/2);
                    Zsize4 = floor(Zsize4/2);  
                    [Xsize5,Ysize5,Zsize5] = size(E5{Axis.x}.array());
                    Xsize5 = floor(Xsize5/2);
                    Ysize5 = floor(Ysize5/2);
                    Zsize5 = floor(Zsize5/2);  
                    %}
%{
                    Edir(1, 1, ci, wi, hi, si, lami) = abs(E1{Axis.x}.array(Xsize1,Ysize1,Zsize1))^2;
                    Edir(1, 2, ci, wi, hi, si, lami) = abs(E1{Axis.y}.array(Xsize1,Ysize1,Zsize1))^2;
                    Edir(1, 3, ci, wi, hi, si, lami) = abs(E1{Axis.z}.array(Xsize1,Ysize1,Zsize1))^2;
                    
                    Edir(2, 1, ci, wi, hi, si, lami) = abs(E2{Axis.x}.array(Xsize2,Ysize2,Zsize2))^2;
                    Edir(2, 2, ci, wi, hi, si, lami) = abs(E2{Axis.y}.array(Xsize2,Ysize2,Zsize2))^2;
                    Edir(2, 3, ci, wi, hi, si, lami) = abs(E2{Axis.z}.array(Xsize2,Ysize2,Zsize2))^2;
                    Edir(3, 1, ci, wi, hi, si, lami) = abs(E3{Axis.x}.array(Xsize3,Ysize3,Zsize3))^2;
                    Edir(3, 2, ci, wi, hi, si, lami) = abs(E3{Axis.y}.array(Xsize3,Ysize3,Zsize3))^2;
                    Edir(3, 3, ci, wi, hi, si, lami) = abs(E3{Axis.z}.array(Xsize3,Ysize3,Zsize3))^2;
                    Edir(4, 1, ci, wi, hi, si, lami) = abs(E4{Axis.x}.array(Xsize4,Ysize4,Zsize4))^2;
                    Edir(4, 2, ci, wi, hi, si, lami) = abs(E4{Axis.y}.array(Xsize4,Ysize4,Zsize4))^2;
                    Edir(4, 3, ci, wi, hi, si, lami) = abs(E4{Axis.z}.array(Xsize4,Ysize4,Zsize4))^2;                    
                    Edir(5, 1, ci, wi, hi, si, lami) = abs(E5{Axis.x}.array(Xsize5,Ysize5,Zsize5))^2;
                    Edir(5, 2, ci, wi, hi, si, lami) = abs(E5{Axis.y}.array(Xsize5,Ysize5,Zsize5))^2;
                    Edir(5, 3, ci, wi, hi, si, lami) = abs(E5{Axis.z}.array(Xsize5,Ysize5,Zsize5))^2;  
                    disp(wi);
                    
                  %  Estrength(1, ci, wi, hi, si, lami) = ...
                   %                   abs(E1{Axis.x}.array(Xsize1,Ysize1,Zsize1))^2 ...
                    %                + abs(E1{Axis.y}.array(Xsize1,Ysize1,Zsize1))^2 ...
                     %               + abs(E1{Axis.z}.array(Xsize1,Ysize1,Zsize1))^2;
   %                 disp(Estrength(1, ci, wi, hi, si, lami));
                    disp(ci);
                    disp(wi);
                    disp(hi); 
                    disp(si);
                    disp(lami);
                                
                    Estrength(2, ci, wi, hi, si, lami) = ...
                                      abs(E2{Axis.x}.array(Xsize2,Ysize2,Zsize2))^2 ...
                                    + abs(E2{Axis.y}.array(Xsize2,Ysize2,Zsize2))^2 ...
                                    + abs(E2{Axis.z}.array(Xsize2,Ysize2,Zsize2))^2;
                    Estrength(3, ci, wi, hi, si, lami) = ...
                                      abs(E3{Axis.x}.array(Xsize3,Ysize3,Zsize3))^2 ...
                                    + abs(E3{Axis.y}.array(Xsize3,Ysize3,Zsize3))^2 ...
                                    + abs(E3{Axis.z}.array(Xsize3,Ysize3,Zsize3))^2;
                    Estrength(4, ci, wi, hi, si, lami) = ...
                                      abs(E4{Axis.x}.array(Xsize4,Ysize4,Zsize4))^2 ...
                                    + abs(E4{Axis.y}.array(Xsize4,Ysize4,Zsize4))^2 ...
                                    + abs(E4{Axis.z}.array(Xsize4,Ysize4,Zsize4))^2;
                    Estrength(5, ci, wi, hi, si, lami) = ...
                                      abs(E5{Axis.x}.array(Xsize5,Ysize5,Zsize5))^2 ...
                                    + abs(E5{Axis.y}.array(Xsize5,Ysize5,Zsize5))^2 ...
                                    + abs(E5{Axis.z}.array(Xsize5,Ysize5,Zsize5))^2;    
                           %}
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
