function [Ls,Rs,thetas] = generate_geometry(obj,L)

    if (obj.verbose)
        display('    generating the geometry ...');
    end
    
    spacings = obj.M*obj.lambda*obj.beta;               % spacings required between adjacent waveguides (m)

    %%
    Nout = ceil(L/spacings);                            % calculate number of needed output waveguides
    Nsplits  = ceil(log(Nout)/log(obj.N_per_split));    % calculate number of splits needed (round up)
    Nout = obj.N_per_split^Nsplits;                     % recalculate new number of output waveguides after rounding
    L_new = Nout*spacings;                              % recalculate new interaction length after rounding

    hs = zeros(Nsplits,1);                              % vertical climbs needed to be accomplished by each bend section
    Rs = zeros(Nsplits,1);                              % optimal radii of curvature calculated (m)
    thetas = zeros(Nsplits,1);                          % optimal angles of bends (rad)
    ds = zeros(Nsplits,1);                              % horizontal distance traversed (m)
    Ls = zeros(Nsplits,1);                              % total length of bent path (m)

    % inline functions to determine angle, path difference for each bend
    % section
    theta_calc = @(h,R) (h < 2*R).*acos(1-h./2./R) + (h >= 2*R).*pi/2;
    d_calc     = @(h,R) 2.*R.*sin(theta_calc(h,R));
    L_calc     = @(h,R) (h < 2*R).*2.*R.*theta_calc(h,R) + (h >= 2*R).*(pi.*R+h-2.*R);
    diff_calc  = @(h,R)  L_calc(h,R)-d_calc(h,R);

    % radius values to try
    NR = 1000000;
    R_list = 10.^(linspace(-6,10,NR));

    %for each bend section, calculate R needed
    for i = (1:Nsplits)    
        h = L_new/obj.N_per_split^i;                                % required vertical climb
        hs(i)      = h;                                             % add to array
        dl_desired = h/obj.beta/obj.n_g;                            % path length difference required for GV matching   
        [~,I] = max(1./abs(diff_calc(h,R_list)-dl_desired));        % solve the path differeces = required for R
        R_opt = R_list(I);                                          % get the R value
        Rs(i) = R_opt;                                              % add to array
        thetas(i) = theta_calc(h,R_opt);
        ds(i) = d_calc(h,R_opt);
        Ls(i) = L_calc(h,R_opt);
    end
    
end