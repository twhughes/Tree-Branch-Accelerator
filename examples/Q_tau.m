
NQ = 50;
Ntau = 51;
NE0 = 20;

Q_list   = 10.^linspace(0,5,NQ);
tau_list = 10.^linspace(-14,-9,Ntau);

factors = zeros(NQ,Ntau);
E0s = zeros(NQ,Ntau);
Gs = zeros(NQ,Ntau);
E_gains = zeros(NQ,Ntau);

s.verbose = false;
     upd = textprogressbar(NQ*Ntau, 'barlength', 20, ...
                               'updatestep', 200, ...
                               'startmsg', '        ',...
                               'endmsg', ' Done! ', ...
                               'showbar', true, ...
                               'showremtime', true, ...
                               'showactualnum', true, ...
                               'barsymbol', '=', ...
                               'emptybarsymbol', ' '); 
s.solve_geometry;

iter = 0;

for qi = (1:NQ)
    s.Q = Q_list(qi);    
    for ti = (1:Ntau)
        iter = iter + 1;
        upd(iter)
        s.tau = tau_list(ti);        
        s.compute_P_max;        
        factors(qi,ti) = s.limiting_factor_index;
        E0s(qi,ti) = s.E0;
        s.compute_gradient;
        Gs(qi,ti) = s.G;
        E_gains(qi,ti) = s.E_gain;
    end
end


%% 
figure(1);
subplot(1,4,4);
imagesc(log10(tau_list*1e12), log10(Q_list), factors);
title('(1=input damage, 2=acc. damage, 3=SPM, 4=SF)');
xlabel('log10 tau (ps)');
ylabel('log10 Q-factor')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
set(gca,'Ydir','normal')
colormap(pmkmp(256,'IsoAZ'));
colorbar();

%{
subplot(3,1,2)
imagesc(log10(tau_list*1e12), log10(Q_list), log10(E0s));
title('log10 E0 (V/m)');
xlabel('log10 tau (ps)');
ylabel('log10 Q-factor')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
set(gca,'Ydir','normal')
colorbar();
colormap(pmkmp(256,'cubicL'));
%}
subplot(1,4,1)
imagesc(log10(tau_list*1e12), log10(Q_list), (Gs/1e6));
title('gradient (MV/m)');
xlabel('log10 tau (ps)');
ylabel('log10 Q-factor')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
set(gca,'Ydir','normal')
colorbar();
colormap(pmkmp(256,'cubicL'));

subplot(1,4,2)
imagesc(log10(tau_list*1e12), log10(Q_list), (E_gains/1e3));
title('energy gain (keV)');
xlabel('log10 tau (ps)');
ylabel('log10 Q-factor')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
set(gca,'Ydir','normal')
colorbar();
colormap(pmkmp(256,'cubicL'));

subplot(1,4,3)
imagesc(log10(tau_list*1e12), log10(Q_list), (E0s)/1e6);
title('maximum E0 (MV/m)');
xlabel('log10 tau (ps)');
ylabel('log10 Q-factor')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','normal')
set(gca,'FontSize',16,'fontWeight','normal')
set(gca,'Ydir','normal')
colorbar();
colormap(pmkmp(256,'cubicL'));
