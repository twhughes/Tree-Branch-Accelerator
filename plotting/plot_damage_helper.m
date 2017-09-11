function [] = plot_damage_helper(obj)
    NT = 1000;
    taus = 10.^linspace(-15,-10,NT);
    Eds = zeros(NT,1);
    for ti = (1:NT)
       Eds(ti) = obj.Ed(taus(ti)); 
    end
    
    
    plot(taus*1e15,Eds/1e9);
    xlabel('pulse duration (fs)')
    ylabel('damage threshold electric field (GV/m)');
    set(gca, 'Xscale','log');
    title(['material stack = ', obj.material_stack]);
    grid on;
end