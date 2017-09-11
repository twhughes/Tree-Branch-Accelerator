function [R_best] = optimize_pillar_radius(obj)
    
    if (obj.verbose)
        display('    optimizing pillar radius ... ');
        upd = textprogressbar(obj.NR, 'barlength', 20, ...
                             'updatestep', 1, ...
                             'startmsg', '        ',...
                             'endmsg', ' Done!', ...
                             'showbar', true, ...
                             'showremtime', true, ...
                             'showactualnum', true, ...
                             'barsymbol', '=', ...
                             'emptybarsymbol', ' ');        
    end
    
    f0 = obj.c0/obj.lambda;                                 % central frequency (Hz)    
    obj.f0 = f0;
    rs = linspace(1e-7,obj.beta*obj.lambda/2,obj.NR);       % radius list to check (m)
                     
    Gs_radius = zeros(obj.NR,1);                            % gradients to store (E0)
    for ri = (1:obj.NR);
        if (obj.verbose)
            upd(ri);
        end        
        r = rs(ri);                                 
        [g,~,fields,ER,E0] = solve_dual_pillars(obj, r, f0);   % solve for gradient (E0)
        Gs_radius(ri) = g;
    end
    obj.Gs_radius = Gs_radius;
    obj.radius_list = rs;
    [G_best,I_best] = max(abs(Gs_radius));
    R_best = rs(I_best);    
    obj.T0 = G_best;
    
    [~,~,fields,ER] = solve_dual_pillars(obj, R_best, f0);
    obj.DLA_fields = fields;
    obj.DLA_struct = ER;
    
    if (obj.verbose)
        display(['        radius   =  ', num2str(R_best/1e-6),' um.'])
        display(['        gradient =  ', num2str(G_best), ' E0.']);   
    end    
end