
NE = 20;
NC = 21;

E0s = 10.^linspace(8,10,NE);
Cs = linspace(-20,20,NC);

Gs = zeros(NE,NC);
damage_map = zeros(NE,NC);

%s.L0 = 1e-2;
%s.L0 = 1e-3;
%s.L0 = 1e-4;
s.L0 = 1e-5;

s.setup;
upd = textprogressbar(NE*NC);

s.verbose = false;
iter = 0;
for ei = (1:NE)
    s.E0 = E0s(ei);
    for ci = (1:NC);
        iter = iter + 1;
        upd(iter)
        s.chirp = Cs(ci);
        s.construct_initial_pulse;       
        s.run_simulation;
        
        Gs(ei,ci) = s.G;
        damage_map(ei,ci) = s.damaged;
    end
end

%%
figure(1)
imagesc(Cs,linspace(5,10,NE),Gs)
xlabel('chirp parameter')
ylabel('incident peak electric field')
title('log10 Gradient. SOI');
colorbar();

figure(2)
imagesc(Cs,linspace(5,10,NE),damage_map)
xlabel('chirp parameter')
ylabel('incident peak electric field')
title('log10 Gradient. SOI');