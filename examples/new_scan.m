

s = simulation('SiN');
s.M = 4;
s.setup();

s.E0 = s.Ed(s.tau);
s.verbose = false;

NL = 100;
Ls = 10.^linspace(-5,1,NL);

Gs = zeros(NL,1);
Es = zeros(NL,1);

upd = textprogressbar(NL);

for i = (1:NL)
    upd(i);
    s.L = Ls(i);
    s.solve_geometry;
    s.run_simulation;
    Gs(i) = s.G;
    Es(i) = s.G*Ls(i);
end