function [A_t_out] = propagate_dist_NLSE(obj,L_dist)

    %% constants and inputs
    c = obj.c0;                         % speed of light (m/s)
    T0 = obj.tau;                       % laser duration (FWHM) (s)
    TR = obj.TR;                        % raman time constant (s)
    lambda0 = obj.lambda;               % laser wavelengh (m)
    n2 = obj.n2;                        % nonlinear index (m^2/W)
    Aeff = obj.A_eff;                   % effective area (m^2)
    gamma = 2*pi/lambda0*n2/(Aeff);     % nonlinear parameter per unit length and power, https://www.rp-photonics.com/self_phase_modulation.html
    E0 = max(abs(obj.A_t));             % maximum E_field in time domain (V/m)
    P0 = E0^2/(2*obj.Z0)*Aeff;          % peak power (W)
    
    % dispersion and nonlinearity parameters
    beta2 = obj.beta2;                  % GVD parameter (fs^2/mm) (GVD)
    beta2 = beta2*1e-27;                % GVD parameter (s^2/m) (GVD)
    distance = L_dist;                  % distance to propagate (m)
    beta2_n = sign(beta2);              % normal/anomalous dispersion (+- 1)
    N = sqrt(gamma*P0*T0^2/abs(beta2)); % nonlinear parameter (unitless)
    
    % simulation parameters
    nt = obj.NT;                        % number of FFT points
    Tmax = obj.T_max/obj.tau;           % maximum time (in units of pulse duration)
    step_num = 20;                      % No. of z steps
    deltaz = distance/step_num;         % step size in z (m)
    dtau = (Tmax)/nt;                   % step size in time (in units of pulse duration)
    omega = (pi/Tmax) * [(0:nt/2-1) (-nt/2:-1)];    % frequency grid
 
    %% Main propagation loop (email Si Tan if you have questions)
    
    A0 = obj.A_t/E0;                                % normalized initial pulse (time domain)    
    dispersion = exp(0.5i*beta2_n*omega.^2*deltaz); % dispersion phase factor
    hhz = 1i*N^2*deltaz;                            % nonlinear phase factor   
    temp = A0.*exp(abs(A0).^2.*hhz/2);              % note hhz/2
        
    for n=1:step_num
        f_temp = ifft(temp).*dispersion;
        A = fft(f_temp);
        temp = A.*exp((abs(A).^2+...
            1*1i/(2*pi/lambda0*c*T0).*(([diff((abs(A).^2).*A) 0]./dtau)./(A+1e-40))+...
            1*TR/T0*[diff((abs(A).^2)) 0]./dtau).*hhz);
    end
    
    A = temp.*exp(-(abs(A).^2+...
        1*1i/(2*pi/lambda0*c*T0).*(([diff((abs(A).^2).*A) 0]./dtau)./(A+1e-40))+...
        1*TR/T0*[diff((abs(A).^2)) 0]./dtau).*hhz/2);

    A_t_out = A*E0;                                 % output pulse
    
end
