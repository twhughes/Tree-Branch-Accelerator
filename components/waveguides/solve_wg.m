function [n_eff, n_g, beta2, A_eff] = solve_wg(obj, wg_width, wg_height, clad_height, n_core, n_clad_top, n_clad_bottom)
    
    if (obj.verbose)        
        display('    solving waveguide parameters...');
    end
    
    % Convert everything into microns
    lambda_um = obj.lambda*1e6;
    wg_width_um = wg_width*1e6;
    wg_height_um = wg_height*1e6;
    clad_height_um = clad_height*1e6;
    
    % Horizontal dimensions:
    rw = wg_width_um/2;           % Ridge half-width 5.3
    side = 0.2;                   % Space on side

    % Grid size:
    %{
    if strcmp(obj.material_stack,'SOI')
        dx = 0.0125/2;            % grid size (horizontal)
        dy = 0.0125/2;            % grid size (vertical)
    else
        %}
        dx = 0.0125/2;              % grid size (horizontal)
        dy = 0.0125/2;              % grid size (vertical)        
    %{end%}
    
    nmodes = 1;                   % number of modes to compute
    % create epilon mesh
    [~,~,~,~,~,~,eps,~] = waveguidemesh([n_clad_bottom,n_core,n_clad_top],[clad_height_um,wg_height_um,clad_height_um], ...
                                                wg_height_um,rw,side,dx,dy); 

    % solve the mode
    [Ex,n_eff] = svmodes(lambda_um,n_core,nmodes,dx,dy,eps,'000S','EX');

    %% Effective Index and N2 Calculations
    
    % Total A_eff
    numerator = (sum(sum(abs(Ex).^2))*dx*dy)^2;
    denominator = sum(sum(abs(Ex).^4))*dx*dy;    % https://www.rp-photonics.com/effective_mode_area.html
    A_eff = 2*numerator/denominator*1e-12;       % A_effective (m) (factor of 2 since we are doing a half simulation)
    
    ns = sqrt(eps); % refractive index map
    obj.A_eff_list = [];
    obj.n2_list = [];
    
    % Top cladding A_eff        
    numerator = (sum(sum(abs(Ex.*(ns==obj.n_clad_top)).^2))*dx*dy)^2;
    denominator = sum(sum(abs(Ex.*(ns==obj.n_clad_top)).^4))*dx*dy;
    obj.A_eff_list = [obj.A_eff_list 2*numerator/denominator*1e-12];       % A_effective (m) (factor of 2 since we are doing a half simulation)
        
    % Core cladding A_eff        
    numerator = (sum(sum(abs(Ex.*(ns==obj.n_core)).^2))*dx*dy)^2;
    denominator = sum(sum(abs(Ex.*(ns==obj.n_clad_top)).^4))*dx*dy;
    obj.A_eff_list = [obj.A_eff_list 2*numerator/denominator*1e-12];       % A_effective (m) (factor of 2 since we are doing a half simulation)
          
    % Bottom cladding A_eff    
    numerator = (sum(sum(abs(Ex.*(ns==obj.n_clad_bottom)).^2))*dx*dy)^2;
    denominator = sum(sum(abs(Ex.*(ns==obj.n_clad_top)).^4))*dx*dy;
    obj.A_eff_list = [obj.A_eff_list 2*numerator/denominator*1e-12];       % A_effective (m) (factor of 2 since we are doing a half simulation)            
    
    % normalize so that sum of A_eff_list = A_eff
    obj.A_eff_list = obj.A_eff_list*A_eff/sum(obj.A_eff_list);         % add to A_eff total
    obj.wg_mode_profile = [fliplr(transpose(Ex)), Ex'];

    % lists of refractive index and corresponding n2 values.  For
    % simplicity in loading
    n_list_total  = [obj.nAir obj.nSi  obj.nSiO2  obj.nSi3N4];
    n2_list_total = [0        obj.n2Si obj.n2SiO2 obj.n2Si3N4];
    obj.n_list = [obj.n_clad_top obj.n_core obj.n_clad_bottom];
    
    % load n2 values
    for m = (1:length(obj.n_list))        
        obj.n2_list = [obj.n2_list n2_list_total(n_list_total == obj.n_list(m))];
    end
    
    %% GVD
    Nl = 5;                                                     % number of wavelengths
    lambda_um_list = linspace(lambda_um-0.1,lambda_um+0.1,Nl);  % wavelength list
    
    % get dispersion curves
    
    % Si dispersion (ref: http://aip.scitation.org/doi/10.1063/1.555624)    
    if strcmp(obj.material_stack,'SOI')
        n_data = load('Si.txt');         
        wl = n_data(:,1);
        n = n_data(:,2);
        fittedmodel = fit(wl, n, 'fourier8');
        n_Si = feval(fittedmodel,lambda_um_list)';
    end
    
    % SiO2 dispersion (https://refractiveindex.info/?shelf=main&book=SiO2&page=Malitson)
    n_SiO2 = sqrt(1+0.6961663*lambda_um_list.^2./(lambda_um_list.^2-0.0684043^2)+...
        0.4079426*lambda_um_list.^2./(lambda_um_list.^2-0.1162414^2)+...
        0.8974794*lambda_um_list.^2./(lambda_um_list.^2-9.896161^2)); % SiO2

    % Si3N4 dispersion (https://refractiveindex.info/?shelf=main&book=Si3N4&page=Luke)
    n_Si3N4 = sqrt(1+3.0249*lambda_um_list.^2./(lambda_um_list.^2-0.1353406^2)+...
        40314*lambda_um_list.^2./(lambda_um_list.^2-1239.842^2)); % Si3N4
    
    if strcmp(obj.material_stack,'SOI')
        n_clad_top_list = ones(1,Nl);
        n_core_list = n_Si;
        n_clad_bottom_list = n_SiO2;
    else
        n_clad_top_list = n_SiO2;
        n_core_list = n_Si3N4;
        n_clad_bottom_list = n_SiO2;
    end
    
    % status bar
    if (obj.verbose)
        upd = textprogressbar(Nl, 'barlength', 20, ...
                               'updatestep', 1, ...
                               'startmsg', '        ',...
                               'endmsg', ' Done! ', ...
                               'showbar', true, ...
                               'showremtime', true, ...
                               'showactualnum', true, ...
                               'barsymbol', '=', ...
                               'emptybarsymbol', ' '); 
    end
    
    % scan wavelength
    n_list_total = zeros(Nl,1);       
    for i = 1:Nl        
        if (obj.verbose) 
            upd(i);        
        end
        [~,~,~,~,~,~,eps,~] = waveguidemesh([n_clad_bottom_list(i),n_core_list(i),n_clad_top_list(i)],[clad_height_um,wg_height_um,clad_height_um], ...
                                                wg_height_um,rw,side,dx,dy);

        [~,n_eff] = svmodes(lambda_um_list(i),n_core_list(i),nmodes,dx,dy,eps,'000S','EX');
        n_list_total(i) = n_eff; 
        
    end

    % polynomial fits and calculate parameters
    
    p = polyfit(lambda_um_list,n_list_total',4);
    
    % http://ticc.mines.edu/csm/wiki/images/f/f0/UFO05-Dispersion.pdf
    beta2 = lambda_um.^3/(2*pi*obj.c0^2).*(12*p(1)*lambda_um.^2 + 6*p(2)*lambda_um + 2*p(3))*1e21;       % fs^2/mm    
    GV = obj.c0./(n_eff-lambda_um.*(4*lambda_um.^3*p(1)+3*lambda_um.^2*p(2)+2*lambda_um*p(3)+p(4)));    
    n_g = obj.c0/GV; 
    
    if (obj.verbose)
        display(['        n_eff = ', num2str(n_eff)]);
        display(['        n_g   = ', num2str(n_g)]);
        display(['        beta2 = ', num2str(beta2), ' (fs^2/mm)']);
        display(['        A_eff = ', num2str(A_eff*1e12), ' (um^2)']);
    end
end