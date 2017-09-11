function [g,U,fields,ER,E0] = solve_dual_pillars(obj,r,f)
        
    lambda_p = obj.c0/f;    
    eps_pil = obj.nPillars^2;
    eps_wg  = obj.n_eff^2;
    eps_core = obj.n_core^2;
    
    % resolution
    dl =   obj.lambda/obj.grids_per_lam;
    npml = obj.npml;
    
    % problem dimensions
    spc_src = lambda_p/1;     % spacing between PML and source
    spc_l1 = lambda_p/1;      % spacing between source and waveguide end
    spc_l2 = lambda_p/2;      % spacing between waveguide end and pillar

    % calculate grid points spacings
    pts_src    = round(spc_src/dl);
    pts_l1     = round(spc_l1/dl);
    pts_l2     = round(spc_l2/dl);
    pts_radius = round(r/dl);
    pts_pil = 2*pts_radius;
    pts_gap = round(obj.gap_spc/dl);

    % calculate grid dimensions
    Nx = round(obj.M*obj.beta*lambda_p/dl);
    Ny = npml+pts_src+pts_l1+pts_l2+pts_pil+pts_gap+pts_pil+pts_l2+pts_l1+pts_src+npml; 
    nx = round(Nx/2);
    ny = round(Ny/2);

    % initialize FDFD arrays and parameters
    BC = [-1,-1];
    Pol = 'Hz';
    NPML = [0,0,npml,npml];
    RES = [dl,dl];
    ER  = ones(Nx,Ny);
    MuR = ones(Nx,Ny);
    b  = zeros(Nx,Ny);    
    
    % create source
    b(nx,npml + pts_src) = 1;
    b(nx,npml + pts_src + pts_l1 + pts_l2 + 2*pts_radius + pts_gap + 2*pts_radius + pts_l2 + pts_l1) = -1;
    
    % draw waveguide
    x_wg_up   = round(Nx/2 + obj.wg_width/2/dl);
    x_wg_down = round(Nx/2 - obj.wg_width/2/dl);
    ER(x_wg_down:x_wg_up,:) = eps_wg;
    ER_wg_only = ER;

    % draw gap
    ER(:,npml+pts_src+pts_l1:npml+pts_src+pts_l1+pts_l2+pts_pil+pts_gap+pts_pil+pts_l2,:) = 1;

    % draw pillars
    pil_left_y  = npml + pts_src + pts_l1 +  + pts_l2 + pts_radius;
    pil_right_y = npml + pts_src + pts_l1 +  + pts_l2 + 2*pts_radius + pts_gap + pts_radius;
    
    for i = (1:obj.M)
        pil_pos_x_i = round((2*i-1)*obj.beta*lambda_p/2/dl);    
        for j = (1:Nx)
            x_pos = (j);
            for k = (1:Ny)
                y_pos = (k);            
                if ( sqrt((x_pos-pil_pos_x_i)^2 + (y_pos-pil_left_y)^2) < pts_radius )
                    ER(j,k) = eps_pil;
                end            
                if ( sqrt((x_pos-pil_pos_x_i)^2 + (y_pos-pil_right_y)^2) < pts_radius )
                    ER(j,k) = eps_pil;                
                end
            end
        end
    end
    
    % normalization
    [fields,~] = FDFD(ER_wg_only,MuR,RES,NPML,BC,lambda_p,Pol,b);
    Ex = fields.Ex;
    Ey = fields.Ey;
    E_abs = sqrt(abs(Ex).^2 + abs(Ey).^2);
    E0 = E_abs(round(nx),round(ny));
    
    % gradient calculation    
    [fields,~] = FDFD(ER,MuR,RES,NPML,BC,lambda_p,Pol,b);
    Ex = fields.Ex;
    Ey = fields.Ey;
    I  = abs(Ex).^2 + abs(Ey).^2;
    I(:,1:npml+pts_src+pts_l1+pts_l2+2*pts_radius) = 0;
    I(:,npml+pts_src+pts_l1+pts_l2+2*pts_radius+pts_gap:end) = 0;
    
    
    gap_fields = Ex(:,ny);    
    g = (sum(exp(1i*2*pi*(1:Nx)*obj.M/Nx).*transpose(gap_fields)))/E0/Nx;
    U = sum(sum(ER.*I));
    %U = I(nx,ny);
end