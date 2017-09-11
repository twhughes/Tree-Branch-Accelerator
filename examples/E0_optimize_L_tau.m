s = simulation('SiN');
%s.Q = 16.4956;
%s.T0 = 0.24981;
%s.R_pillars = 1.6122e-6;
s.M = 4;
s.setup;

%%
NE0 = 45;
Ntau = 41;
NL = 40;

E0_list = 10.^linspace(5,10,NE0);
tau_list = 10.^linspace(-15,-12,Ntau);
L_list = 10.^linspace(-5,0,NL);


Gs = zeros(Ntau,NL);
Es = zeros(Ntau,NL);
E0s_max = zeros(Ntau,NL);
Gs_dam = zeros(Ntau,NL);
Es_dam = zeros(Ntau,NL);
E0s_max_dam = zeros(Ntau,NL);

upd = textprogressbar(NE0*Ntau*NL);

iter = 0;
s.verbose = false;
for li = (1:NL)
    s.L = L_list(li);
    s.solve_geometry;
    for ti = (1:Ntau)
        s.tau = tau_list(ti);
        G_local = zeros(NE0,1);
        G_local_dam = zeros(NE0,1);
        
        for ei = (1:NE0)
            iter = iter + 1;
            upd(iter);            
            s.E0 = E0_list(ei);
            s.construct_initial_pulse;
            
            s.run_simulation;
            G_local(ei) = s.G;         
            G_local_dam(ei) = s.G*(1-s.damaged);               
        end
        [G,I] = max(G_local);  
        Gs(ti,li) = G;
        Es(ti,li) = G*s.L;
        E0s_max(ti,li) = E0_list(I);

        [G,I] = max(G_local_dam);  
        Gs_dam(ti,li) = G;
        Es_dam(ti,li) = G*s.L;
        E0s_max_dam(ti,li) = E0_list(I);        
    end
end

%%
figure(1);
imagesc(log10(tau_list),log10(L_list),log10(Gs)'); colorbar;
title('log10 Gradient V/m (without damage)')
xlabel('log10 tau (s)');
ylabel('log10 L (m)');
set(gca,'Ydir','normal')
figure(2);
imagesc(log10(tau_list),log10(L_list),log10(Es)'); colorbar;
title('log10 Energy gain eV (without damage)')
xlabel('log10 tau (s)');
ylabel('log10 L (m)');
set(gca,'Ydir','normal')
figure(3);
imagesc(log10(tau_list),log10(L_list),log10(E0s_max)'); colorbar;
title('log10 E0 V/m (without damage)')
xlabel('log10 tau (s)');
ylabel('log10 L (m)');
set(gca,'Ydir','normal')


figure(4);
imagesc(log10(tau_list),log10(L_list),log10(Gs_dam)'); colorbar;
title('log10 Gradient V/m (with damage)')
xlabel('log10 tau (s)');
ylabel('log10 L (m)');
set(gca,'Ydir','normal')
figure(5);
imagesc(log10(tau_list),log10(L_list),log10(Es_dam)'); colorbar;
title('log10 Energy gain eV (with damage)')
xlabel('log10 tau (s)');
ylabel('log10 L (m)');
set(gca,'Ydir','normal')
figure(6);
imagesc(log10(tau_list),log10(L_list),log10(E0s_max_dam)'); colorbar;
title('log10 E0 V/m (with damage)')
xlabel('log10 tau (s)');
ylabel('log10 L (m)');
set(gca,'Ydir','normal')