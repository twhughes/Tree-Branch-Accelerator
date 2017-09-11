function [] = plot_DLA_helper(obj)

    [Nx,Ny] = size(obj.DLA_struct);
    dl = obj.lambda/obj.grids_per_lam;
    xrange = linspace(0,dl*Nx*1e6,Nx);
    yrange = linspace(0,dl*Ny*1e6,Ny);            
    figure; clf; hold all;
    imagesc(yrange,xrange,obj.DLA_struct);        
    axis equal tight; colorbar();
    xlabel('y position (um)');
    ylabel('x position (um)');
    title('device relative permittivity');

    figure; clf; hold all;
    Hz = obj.DLA_fields.Hz;
    imagesc(yrange,xrange,abs(Hz));        
    axis equal tight; colorbar();
    xlabel('y position (um)');
    ylabel('x position (um)');
    title('abs(Hz)');            

    figure; clf; hold all;
    Ex = obj.DLA_fields.Ex;
    imagesc(yrange,xrange,real(Ex));        
    axis equal tight; colorbar();
    xlabel('y position (um)');
    ylabel('x position (um)');
    title('real(Ex)');            


    figure; clf; hold all;
    Ey = obj.DLA_fields.Ex;
    imagesc(yrange,xrange,real(Ey));        
    axis equal tight; colorbar();
    xlabel('y position (um)');
    ylabel('x position (um)');
    title('real(Ey)');            

    figure; clf; hold all;
    plot(xrange,real(Ex(:,round(Ny/2))));
    xlabel('position along x gap');
    ylabel('real(Ex)');
    title('parallel gap fields');

    figure; clf; hold all;
    plot(yrange,abs(Ex(round(Nx/2),:)));
    xlabel('position perpendicular to gap');
    ylabel('abs(Ex)');
    title('perpendicular gap fields');      

    if (~isempty(obj.Ts_raw) && ~isempty(obj.fs_scan))
        figure; clf; hold all;
        plot(obj.fs_scan/1e12,abs(obj.Ts_raw));
        xlabel('frequency (THz)');
        ylabel('G(E0)');
        title('raw FDFD scan'); 
    end            
    if (~isempty(obj.Ts_fit) && ~isempty(obj.fs_centered))
        hold all;
        plot(obj.fs_centered/1e12,obj.Ts_fit);
        xlabel('frequency (THz)');
        ylabel('G(E0)');
        title('raw FDFD scan'); 
        legend({'raw scan', 'Lorenzian fit'});
    end
    xlim([obj.fs_scan(1)/1e12, obj.fs_scan(end)/1e12])
    
    if (~isempty(obj.Gs_radius) && ~isempty(obj.radius_list))
        
        figure; clf;
        plot(obj.radius_list*1e6,abs(obj.Gs_radius)');
        xlabel('pillar radius (um)')
        ylabel('acc. gradient E_0');
    end
end