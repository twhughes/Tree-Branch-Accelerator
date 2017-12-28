classdef simulation < dynamicprops
    % Class used for simulating tree branch laser-couping structure for dielectric laser accelerators.      
    % 
    %     In this documentation, the 'obj' in 'obj.x' refers to the simulation object that is created by obj = simulation.    
    %     The 'x' refers to either the parameter (below) or the method (further below)
    % 
    %     Here is an example of a use case:
    %         s = simulation;
    %         s.tau = 250e-9;
    %         s.L0 = 1e-3;
    %         s.E0 = 1e9;
    %         
    %         s.solve_waveguide_parameters;
    %         s.load_geometry;
    %
    %         s.compute_P_max;
    %         s.compute_gradient
    %
    %         display(['gradient of ', num2str(s.G*1e6), ' MV/m']);
    
    properties (Constant)        
        % constants        
        c0 = 299792458;                     % speed of light (m/s)
        e0 = 8.854187817e-12;               % vacuum permittivity (F/m)
        Z0 = 376.730313416;                 % impedance of free space (Ohms)
        nAir = 1;                           % refractive index of air   at 2um
        nSi  = 3.451;                       % refractive index of Si    at 2um
        nSiO2 = 1.4381;                     % refractive index of SiO2  at 2um
        nSi3N4 = 1.9834;                    % refractive index of Si3N4 at 2um
        n2Si    = 2e-18;                    % nonlinear index (m^2/W)
        n2SiO2  = 2.6e-20;                  % nonlinear index (m^2/W)
        n2Si3N4 = 2.4e-19;                  % nonlinear index (m^2/W)
    end
    
    properties     
        % Input pulse default parameters
        lambda = 2e-6;                      % wavelength (m)
        w0;
        tau = 250e-15;                      % pulse duration (s)
        E0 = 1e9;                           % input field strength (V/m)
        chirp = 0;                          % input pulse chirp (?)
        tbp = 0.44;                         % time bandwidth product
        rep_rate = 1e6;                     % repitition rate (Hz)
        NT = 2^15;                          % number of pulse time points 
        T_max = 40e-12;                     % maximum time span to look at (s)
        num_taus_sample = 500;              % number of pulse durations to sample for input pulse      
        % Input coupler default parameters
        coupler_efficiency = 0.5;           % input coupler maximum power efficiency
        split_efficiency = 0.95;            % splitting power efficiency
        bend_efficiency = 0.95;             % bending power efficiency
        bend_efficiency_list;               % list of bending power efficiency per bend
        % Waveguide parameters
        material_stack = 'SOI';             % if 'SOI', silicon on insulator.  if 'SiN', Si3N4 core with SiO2 cladding
        n_core;                             % refractive index of waveguide core
        n_clad_top;                         % refractive index of waveguide top cladding
        n_clad_bottom;                      % refractive index of waveguide bottom cladding        
        wg_height = 220e-9;                 % height of waveuide core (m)
        wg_width = 760e-9;                  % width of waveguide core (m)
        clad_height = 1e-6;                 % height of cladding
        n_list;                             % list of refractive index of materials
        n2_list;                            % n2 of waveguide materials (m^2/W)
        A_eff_list;                         % effective area (defined by power) of each material in n2_list.
        n_eff;                              % effective index of waveguide (calculated)
        n_g;                                % group index of waveguide (calculated)
        beta2;                              % group velocity dispersion (fs^2 / mm)
        A_eff;                              % waveguide mode effective area (m^2)
        wg_mode_profile;                    % Ex mode profile (V/m) for plotting
        % Geometry parameters
        beta = 1;                           % electron speed / speed of light
        M = 3;                              % periods of DLA powered by one waveguide
        L = 192e-6;                           % interaction length (m)
        L0 = 1e-5;                          % space between input coupler and first split (m)
        Ls;                                 % list of bend lengths
        Rs;                                 % list of bend radii of curvature
        thetas;                             % list of bend angles
        N_per_split = 2;                    % number of outputs per split (only 2 is functional) 
        % Damage thresholds and nonlinearities
        phi_SPM_limit = 2*pi;                 % maximum allowable SPM phase
        Ed;                                 % anonymous function giving peak field damage threshold (V/m) vs tau (s)        
        P0_SPM;                             % max input power (W) before SPM
        P0_SF;                              % max input power (W) before self-focusing
        E0_SPM;                             % max input field (V/m) before SPM
        E0_SF;                              % max input field (V/m) before self-focusing
        E0_input_dam;                       % input field (V/m) causing damage at input;
        E0_acc_dam;                         % input field (V/m) causing damage at accelerator;
        limiting_factor;                    % string describing the limiting factor
        limiting_factor_index;              % identifying index for factor (1-> input 2-> acc 3-> SPM 4-> SF)
        % DLA parameters
        Q = 120;                            % DLA quality factor
        T_at_Q1 = 0.0357;                   % gradient at f = f0 and Q=1. For Q > 1, G(f0) = T_at_Q1*sqrt(Q)
        ws;
        g_w;
        g_w_ZX;      
        E0_w;
        G_w;
        G_t0;
        enhancement = 2;                    % E_max/E0 (in material) (computed in obj.solve_DLA_parameters)
        fA = 0.5;                           % G/E_material (acceleration factor)
        G;                                  % final gradient (V/m) (note, can be > 0 if structure is damaged, need to check separately with obj.damaged)        
        E_gain;                             % total energy gain (eV).  obj.G * obj.L
        % simulation parameters
        verbose = true;                     % whether to display progress messages.  Keep on for debugging.  Turn off for doing scans.
    end
    
    methods
        function obj = simulation
            % optional argument to constructor 'wg_stack' can be 'SOI' or
            % 'SiN'.  Sets the material stack. If none given, SOI is used
            % as default.  Also loads default wg_geometry and materials for
            % each material stack
            if (obj.verbose)
                display('starting simulation...');
            end
            obj.w0 = 2*pi*obj.c0/obj.lambda;
            
        end
    %---- SHORTCUT METHODS ----%
        function obj = setup(obj)
            % solves the waveguides, geometry, and DLA, loads damage data, constructs initial pulse.
            obj.solve_waveguide_parameters;          
            obj.solve_geometry;
        end
        function obj = run_simulation(obj)
            % performs input coupling, pulse propagation, gradient computation.  
            obj.compute_P_max;
            obj.compute_gradient;
        end
        function obj = plot_all(obj)
            % runs all of the plotting methods at once
            close all;
            obj.plot_waveguide;
            obj.plot_geometry;
            obj.plot_DLA;
            obj.plot_pulse;
        end
    %---- STRUCTURE DEFINITION METHODS ----%

        function obj = solve_waveguide_parameters(obj)
            % does error checking on wg geometry, then solves for n_eff, n_g, beta2, and A_eff
            if (isempty(obj.material_stack)  || isempty(obj.wg_height) || isempty(obj.wg_width) || isempty(obj.clad_height))
                error('must specify material_stack, wg_height, wg_width, clad_height')
            end
            if strcmp(obj.material_stack,'SOI')    
                obj.n_core = obj.nSi;
                obj.n_clad_top = obj.nSiO2;
                obj.n_clad_bottom = obj.nSiO2;
            elseif strcmp(obj.material_stack,'SiN')
                obj.n_core = obj.nSi3N4;
                obj.n_clad_top = obj.nSiO2;
                obj.n_clad_bottom = obj.nSiO2;
            else
                error('incorrect material stack, must be either SOI or SiN');   
            end
            
            % mode solver
            [n_eff_calc, n_g_calc, beta2_calc, A_eff_calc] = solve_wg(obj, obj.wg_width, obj.wg_height, obj.clad_height, obj.n_core, obj.n_clad_top, obj.n_clad_bottom);
            obj.n_eff = n_eff_calc;
            obj.n_g = n_g_calc;
            obj.beta2 = beta2_calc;
            obj.A_eff = A_eff_calc;
            
            obj.load_damage;            
        end
        
        function obj = solve_geometry(obj)
        % solve for the lengths, radii of curvature, and bend angles    
            if (2*obj.wg_width > obj.beta*obj.lambda*obj.M)
                error('waveguide is too large to simulate with given lambda, M, beta.  Make M larger and try again.');               
            end        
            [Ls_calc,Rs_calc,thetas_calc] = generate_geometry(obj,obj.L);
            obj.Ls     = Ls_calc;
            obj.Rs     = Rs_calc;
            obj.thetas = thetas_calc;
            %{
            SOI_Rs = open('~/Documents/Fan/ACHIP/DLA_Theory_Paper/Simple_Code_DLA_Laser_Couple/data/SOI_Rs.mat');
            SOI_Ts = open('~/Documents/Fan/ACHIP/DLA_Theory_Paper/Simple_Code_DLA_Laser_Couple/data/SOI_Ts.mat');
            SiN_Rs = open('~/Documents/Fan/ACHIP/DLA_Theory_Paper/Simple_Code_DLA_Laser_Couple/data/SiN_Rs.mat');
            SiN_Ts = open('~/Documents/Fan/ACHIP/DLA_Theory_Paper/Simple_Code_DLA_Laser_Couple/data/SiN_Ts.mat');
            
            SOI_Rs = SOI_Rs.SOI_Rs;
            SOI_Ts = SOI_Ts.SOI_Ts;
            SiN_Rs = SiN_Rs.SiN_Rs;
            SiN_Ts = SiN_Ts.SiN_Ts;
                        
            for i = (1:length(obj.Rs));
                if strcmp(obj.material_stack,'SOI')    
                    obj.bend_efficiency_list(i) = interp1(log(SOI_Rs*1e-9),SOI_Ts,log(obj.Rs(i)));
                elseif strcmp(obj.material_stack,'SiN')
                    obj.bend_efficiency_list(i) = interp1(log(SiN_Rs*1e-9),SiN_Ts,log(obj.Rs(i)));                    
                else
                    error('incorrect material stack, must be either SOI or SiN');   
                end 
            end
%            bend_loss_list_dB = transpose(10*log(ones(size(obj.bend_efficiency_list))-obj.bend_efficiency_list)/log(10));               
%            dist_list = obj.Rs*pi/2;
%            bend_loss_dB_per_m = bend_loss_list_dB./dist_list;
%            bend_loss_partial_bend = bend_loss_dB_per_m.*obj.Rs.*obj.thetas;                
%            obj.bend_efficiency_list = ones(size(obj.bend_efficiency_list))' - 10.^(bend_loss_partial_bend/10);                
            %}
        end
        function obj = load_damage(obj)
            % load damage data into function obj.Ed(tau) using the core of the given waveguide materials
            Ed_anon_fn = get_LIDT(obj);                                     % get damage E as a fn. of tau
            obj.Ed = Ed_anon_fn;                                            % save to simulation
        end
        
        
        function obj = compute_P_max(obj)
            % computes the fields corresponding to the limiting constraints
            if (obj.verbose)
                display('    computing maximum input powers...');
            end            
            if (isempty(obj.Ls))
                error('must solve_geometry first');
            end
            if (isempty(obj.n2_list) || isempty(obj.A_eff_list))
                error('must solve_waveguide_parameters first');
            end     
            
            if (isempty(obj.Q) || isempty(obj.enhancement))
                error('must specify Q and enhancement first');
            end   
            
            % precalculate commonly occuring quantities
            n2_eff = sum(obj.n2_list .* obj.A_eff_list)/obj.A_eff;
            L_sum  = sum((obj.bend_efficiency*obj.split_efficiency).^(0:length(obj.Ls)-1)' .*obj.Ls ./ 2.^(0:length(obj.Ls)-1)');
            n_sum  = sum(obj.A_eff_list .* obj.n_list)/obj.A_eff;
            
            % calculate input power to trigger nonlinearities
            obj.P0_SPM = obj.lambda*obj.A_eff/L_sum/n2_eff/2/pi*obj.phi_SPM_limit/obj.coupler_efficiency;
            obj.P0_SF  = obj.lambda*obj.A_eff*obj.n_g/obj.tau/obj.c0/n2_eff*sqrt(4*log(2)/pi)/2*obj.coupler_efficiency;
            
            % convert these powers to fields
            obj.E0_SF = sqrt(2*obj.P0_SF/obj.c0/obj.e0/n_sum/obj.A_eff);
            obj.E0_SPM = sqrt(2*obj.P0_SPM/obj.c0/obj.e0/n_sum/obj.A_eff);
            
            % get input field damage thresholds            
            N_splits = length(obj.Ls)-1;
            
            obj.E0_input_dam = obj.Ed(obj.tau);
            %obj.E0_acc_dam = obj.Ed(obj.tau)*obj.fA/obj.T_at_Q1*(obj.coupler_efficiency*obj.split_efficiency^N_splits*obj.bend_efficiency^N_splits)^(-1/2)/sqrt(obj.Q);
            obj.E0_acc_dam = obj.Ed(obj.tau)/obj.fA*(obj.coupler_efficiency*obj.split_efficiency^N_splits*obj.bend_efficiency^N_splits*2^(-N_splits))^(-1/2)/sqrt(obj.Q);
            
            factors   = {'Input Damage' 'Accelerator Damage' 'SPM' 'SF'};
            E0_limits = [obj.E0_input_dam obj.E0_acc_dam obj.E0_SPM obj.E0_SF];
            [E0_min, I] = min(E0_limits);
            obj.limiting_factor = factors(I);
            obj.limiting_factor_index = I;
            obj.E0 = E0_min;
        end
        
        function obj = compute_gradient(obj)
            % uses Lorentzian fit to approximate frequency response and
            % solve for acceleration gradient and energy gain.
            if (obj.verbose)
                display('    computing gradient...');
            end                    
            N_splits = length(obj.Ls)-1;
            %field in waveguide headed to accelerator
            E0_acc = obj.E0*2^(-N_splits/2)*sqrt(obj.coupler_efficiency)*obj.split_efficiency^(N_splits/2)*obj.bend_efficiency^(N_splits/2);           
            
            
            obj.T_max = obj.tau*obj.num_taus_sample;                                                        % set maximum time
            dt = obj.T_max/obj.NT;                                                          % time step 
            ts = (-obj.NT/2:obj.NT/2-1)*dt;                                                 % set time points
            E0_t = E0_acc*exp(-2*log(2)*(1+1i*obj.chirp).*ts.^2/obj.tau^2);                 % get initial pulse time series                  
            f_sample = 1/dt;                                                                % sampling frequency
            df = f_sample/obj.NT;                                                           % frequency spacing
            fs = -f_sample/2:df:f_sample/2-df;                                              % frequency array
            obj.ws = fs.*2*pi + obj.w0*ones(size(fs));            
            obj.g_w = obj.T_at_Q1*obj.w0/2/sqrt(obj.Q).*( 1i*(obj.ws-obj.w0*ones(size(obj.ws))) + obj.w0/2/obj.Q*ones(size(obj.ws)) ).^(-1);
           % obj.g_w_ZX = obj.T_at_Q1*()*sqrt(obj.Q/4/pi) ...
          %      .*exp(1i*(obj.ws - obj.w0)./obj.w0*2*obj.Q) ...
         %       .*exp(-2*log(2)*(obj.ws-obj.w0).^2*((obj.Q/obj.w0)^2)+(obj.tau/4/log(2))^2);
            obj.E0_w = fftshift(fft(E0_t));
            obj.G_w = obj.g_w.*abs(obj.E0_w);
            obj.G_t0 = (ifft(obj.G_w));
            obj.G = max(abs(obj.G_t0));
            obj.E_gain = obj.G*obj.L;
        end
        
    %---- AUXILARY METHODS (PLOTTING MOSTLY) ---%
        function obj = plot_geometry(obj)
            % draws a figure of the geometry that was solved for
            draw_structure(obj, obj.Ls, obj.Rs, obj.thetas);
        end
        function obj = plot_waveguide(obj)            
            % plots the waveguide mode profile and n_eff dispersion
            if (isempty(obj.material_stack)  || isempty(obj.wg_height) || isempty(obj.wg_width) || isempty(obj.clad_height))
                error('must specify material_stack, wg_height, wg_width, clad_height')
            end
            if strcmp(obj.material_stack,'SOI')    
                obj.n_core = obj.nSi;
                obj.n_clad_top = obj.nAir;
                obj.n_clad_bottom = obj.nSiO2;
            elseif strcmp(obj.material_stack,'SiN')
                obj.n_core = obj.nSi3N4;
                obj.n_clad_top = obj.nSiO2;
                obj.n_clad_bottom = obj.nSiO2;
            else
                error('incorrect material stack, must be either SOI or SiN');   
            end
            draw_waveguides(obj, obj.wg_width, obj.wg_height, obj.clad_height, obj.n_core, obj.n_clad_top, obj.n_clad_bottom);
        end 
        function obj = plot_damage(obj)
            % plot the damage threshold of the structure, as a function of
            % tau
            display('    plotting damage thresholds...');            
            plot_damage_helper(obj);           
        end   
        function [obj, L_final_list] = get_Length_list(obj, L_max)
            % convenience function.  Takes 'L_max' argument as the longest
            % accelerator interaction length you are interested in in a
            % scan.  Returns a list of all of the lengths that need to be
            % considered up to L_max (rounding up for 2^x).  This means you
            % dont have to consider interaction lengths in between these
            % ones in a scan.
            L_section  = obj.M*obj.beta*obj.lambda;         % length of each waveguide-powered section
            N_sections = ceil(L_max/L_section);             % Number of sections needed to cover L_max
            Nsplits = ceil(log(N_sections)/log(2));         % number of 1->2 splits to get N_sections out
            L_final_list = zeros(Nsplits+1,1);              % There are Nsplits+1 different possible lengths up to Nsections 
            for i = (0:Nsplits)
                L_final_list(i+1) = L_section*2^i;          % compute the interaction lengths of each of these ways
            end        
            L_final_list = L_final_list(2:end);
        end
    end
end