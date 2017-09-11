s = simulation;
s.verbose = false;
NE = 20;
NQ = 21;
s.L = 1e-5;
s.Q = 1;
s.T0 = 1;
s.T_at_Q1 = 0.0391;
s.R_pillars = 1e-6;
s.setup;
E0_list = 10.^linspace(4,10,NE);
Q_list = 10.^linspace(0,6,NQ);

upd = textprogressbar(NE*NQ, 'barlength', 20, ...
                     'updatestep', 5, ...
                     'startmsg', '        ',...
                     'endmsg', ' Done! ', ...
                     'showbar', true, ...
                     'showremtime', true, ...
                     'showactualnum', true, ...
                     'barsymbol', '=', ...
                     'emptybarsymbol', ' ');

Gs = zeros(NQ,NE); 
Es = zeros(NQ,NE);
damaged = zeros(NQ,NE);
s.chirp = 0;
iter = 1;
for qi = (1:NQ)
    s.Q = Q_list(qi);
    s.T0 = s.T_at_Q1*sqrt(s.Q);
    s.solve_DLA_parameters;
    for ei = (1:NE)
        upd(iter);
        iter = iter + 1;
        s.E0 = E0_list(ei);
        s.verbose = false;
        s.construct_initial_pulse;
        s.input_couple;
        s.propagate_pulse;
        s.compute_gradient;
        Gs(qi,ei) = s.G;
        Es(qi,ei) = s.G*s.L;
        damaged(qi,ei) = s.damaged;
    end
end

%%
figure(1); clf;
imagesc(log10(E0_list), log10(Q_list),log10(Gs))
title('log10 gradient V/m')
xlabel('log10 E0 (V/m)')
ylabel('log10 Q-factor (m)')

colorbar()
figure(2); clf;
imagesc(log10(E0_list), log10(Q_list), log10(Es))
title('log10 energy gain eV')
xlabel('log10 E0 (V/m)')
ylabel('log10 Q-factor (m)')
colorbar()


colorbar()
figure(3); clf;
imagesc(log10(E0_list), log10(Q_list), damaged);
title('damaged')
xlabel('log10 E0 (V/m)')
ylabel('log10 Q-factor (m)')
colorbar()