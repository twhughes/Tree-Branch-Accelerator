NE = 20;
NL = 1;

E0_list = 10.^linspace(4,10,NE);
L_list = 10.^linspace(-3,-3,NL);

upd = textprogressbar(NE*NL, 'barlength', 20, ...
                     'updatestep', 5, ...
                     'startmsg', '        ',...
                     'endmsg', ' Done! ', ...
                     'showbar', true, ...
                     'showremtime', true, ...
                     'showactualnum', true, ...
                     'barsymbol', '=', ...
                     'emptybarsymbol', ' ');

Gs = zeros(NL,NE); 
Es = zeros(NL,NE);                     
s.chirp = 2;
iter = 1;
for li = (1:NL)
    s.L = L_list(li);
    s.solve_geometry;
    for ei = (1:NE)
        upd(iter);
        iter = iter + 1;
        s.E0 = E0_list(ei);
        s.verbose = false;
        s.construct_initial_pulse;
        s.input_couple;
        s.propagate_pulse;
        s.compute_gradient;
        Gs(li,ei) = s.G;
        Es(li,ei) = s.G*s.L;
    end
end

%%
figure(1); clf;
imagesc(log10(E0_list), log10(L_list),log10(Gs))
title('log10 gradient V/m')
xlabel('log10 E0 (V/m)')
ylabel('log10 Length (m)')

colorbar()
figure(2); clf;
imagesc(log10(E0_list), log10(L_list), log10(Es))
title('log10 energy gain eV')
xlabel('log10 E0 (V/m)')
ylabel('log10 Length (m)')

colorbar()
