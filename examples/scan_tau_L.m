Ntau = 20;
NL = 21;

taus = 10.^linspace(-15,-10,Ntau);
Ls = 10.^linspace(-5,1,NL);

Gs = zeros(Ntau,NL);
E_gains = zeros(Ntau,NL);

upd = textprogressbar(NL*Ntau, 'barlength', 20, ...
                               'updatestep', 1, ...
                               'startmsg', '        ',...
                               'endmsg', ' Done! ', ...
                               'showbar', true, ...
                               'showremtime', true, ...
                               'showactualnum', true, ...
                               'barsymbol', '=', ...
                               'emptybarsymbol', ' ');    

iter = 1;                           
for li = (1:NL)
    L = Ls(li);
    
    for ti = (1:Ntau)
        tau = taus(ti);
            upd(iter)
        E0 = s.Ed(tau);
        s.specify_input_pulse(lambda,E0,tau,tbp,chirp,rep_rate);
        s.construct_initial_pulse();
        s.input_couple();
        s.specify_geometry(L, beta, M, L0);        
        s.solve_geometry();
        s.propagate_pulse();
        s.compute_gradient();
        
        Gs(ti,li) = s.G; s.G;
        E_gains(ti,li) = s.E_gain;
        
        iter = iter + 1;
    end
end

%%
figure(1); clf;
imagesc((Gs))
figure(2); clf;
imagesc((E_gains));