function frequency_scan(obj)
    
    % ensure the frequency grid is set before doing the lorenzian fit
    if (isempty(obj.fs))
        error('please construct_initial_pulse before solving DLA parameters');        
    end
    
    % get releavant frequency data
    f0 = obj.c0/obj.lambda;   
    fs_centered = obj.fs + f0*ones(size(obj.fs));
    obj.fs_centered = fs_centered;
    
    % if no T0 or Q are given, need to do FDFD scan    
    if (isempty(obj.T0) || isempty(obj.Q))
        
        % progress bar if verbose
        if (obj.verbose)
            display('    scanning frequency to compute Q and Ts ...');
            upd = textprogressbar(obj.NF, 'barlength', 20, ...
                             'updatestep', 1, ...
                             'startmsg', '        ',...
                             'endmsg', ' Done! ', ...
                             'showbar', true, ...
                             'showremtime', true, ...
                             'showactualnum', true, ...
                             'barsymbol', '=', ...
                             'emptybarsymbol', ' ');
        end        
        
        nf = round(obj.NF/2);                                   % center frequency index  
        Ts = zeros(obj.NF,1);                                   % gradients
        Us = zeros(obj.NF,1);                                   % energies       
        delta_f = obj.M/obj.tau*2;                              % length-determined characteristic frequency spacing to scan (see paper)
        fs_scan = linspace(f0-delta_f, f0+delta_f, obj.NF);     % frequencies to scan (Hz)
        obj.fs_scan = fs_scan;                                  % load into simulation
        df_scan = (fs_scan(end)-fs_scan(1))/obj.NF;             % get the frequency spacing (Hz)
        % scan frequencies and get the peaks to fit
        for fi = (1:obj.NF)
            if (obj.verbose)
                upd(fi);
            end        
            f = fs_scan(fi);
            [g,U,~,~] = solve_dual_pillars(obj, obj.R_pillars, f);
            Ts(fi) = g;
            Us(fi) = U;
        end            
        obj.Ts_raw = Ts;
        
        % fit the gradient^2 spectrum to lorenzian        
        fs_norm = (1:obj.NF);                               % normalized frequency grid
        T_2 = abs(Ts).^2;                                   % gradient^2 
        options = optimset('lsqcurvefit');                  % least squares curve fit
        options.Display = 'off';                            % no display

        P0 = [T_2(nf) nf 1];                                % starting guesses [P1 P2 P3] for lorenzian of form F(X) = P1/((X-P2)^2 + P3)
        BOUNDS = [0 nf-nf/3 -inf; inf nf+nf/3 inf];         % bounds for the fit parameters
        [T2FIT, fit_params, ~, ~] = lorentzfit(fs_norm,T_2',P0,BOUNDS,'3',options);    % do the lorenzian fit
        P2 = fit_params(2);                                 % get P2 (center of peak)
        P3 = fit_params(3);                                 % get P3
        TFIT = sqrt(T2FIT);                                 % T2FIT is the fit of gradient^2, TFIT is the fit of gradient
        f0_prime = interp1(fs_norm, fs_scan, P2);           % interpolate the normalized grid with the fs_scan to get the peak center in Hz
        FWHM = 2*sqrt(P3);                                  % compute FWHM of the lorenzian (in normalized frequency)
        FWHM_f = FWHM*df_scan;                              % compute FWHM in Hz
        Q = f0_prime/FWHM_f;                                % calculate the resulting Q-factor
        obj.Q = Q;                                          % load into simulation
        obj.f0_prime = f0_prime;                            % load into simulation
        
        % transform fitted tranfer function |T(w)|^2 in Hz frequency xaxis               
        T0 = max(TFIT);                                     % get the maximum (T0)
        obj.T0 = T0;                                        % load into simulation

        P2 = obj.f0_prime;                                  % new P2
        P3 = obj.f0_prime^2/4/obj.Q^2;                      % new P3
        P1 = T0^2*P3;                                       % new P1
        
        Ts_2_fit = P1./((fs_centered - P2*ones(size(fs_centered))).^2 + P3);    % calculate lorenzian in full frequecny grid
        Ts_fit = sqrt(Ts_2_fit);                            % take sqrt() to get T(w)
        obj.Ts_fit = Ts_fit;                                % load into simulation
        
        obj.T_at_Q1 = T0/sqrt(Q);                           % this quantity represents the gradient at f0 for a Q = 1.  For Q > 1, G(f0) = T_at_Q1*sqrt(Q) 
        % NOTE: this allows us to estimate the response for higher Q-factors.
        % Taking into account the higher T0 but also smaller bandwidth.
        
        if (obj.verbose)
            display(['        Q  = ', num2str(Q)  ]);
            display(['        T0 = ', num2str(T0)  ]);            
        end        
        
    else
        
        % In this case, we are given a Q and T0 to work with already
        if (obj.verbose)
            display(['    using specified Q = ', num2str(obj.Q), ' and T0 = ', num2str(obj.T0),' to approximate DLA response']);
        end
        
        % recover the fitted transfer function with Q and T0        
        if isempty(obj.f0_prime)
            obj.f0_prime = f0;
        end
        P2 = obj.f0_prime;
        P3 = obj.f0_prime^2/4/obj.Q^2;
        P1 = obj.T0^2*P3;
        Ts_2 = P1./((fs_centered - P2*ones(size(fs_centered))).^2 + P3);
        obj.Ts_fit = sqrt(Ts_2);
        
    end 

    % compute field enhancement in structure
    [~,~,fields,ER,E0] = solve_dual_pillars(obj, obj.R_pillars, f0);
    E_abs = sqrt(abs(fields.Ex/E0).^2 + abs(fields.Ey/E0).^2);
    E_abs_in_structure = (ER > 1).*E_abs;
    obj.enhancement = max(E_abs_in_structure(:));
end
