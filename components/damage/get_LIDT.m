function [Ed] = get_LIDT(obj)

    % Damage Data
    % "Laser damage threshold measurements of optical materials for direct
    % laser accelerators"
    data_SiO2 = [
        -3,  0.8;
        -2,  1.2;
        -1,  1.8; 
        0,   2.2; 
        1,   5; 
        2,   12; 
        3,   40
    ];  %log10ps -> J/cm^2

    % Ken Soong's damge fluence data (for scaling)
    SiO2_factor = 3.5;
    Si_factor = 0.20;                                        
    Si3N4_factor = 0.65;                                             % needs confirmation

    durations = data_SiO2(:,1);                                      % Pulse durations log10(ps)
    Fs = data_SiO2(:,2);                                             % Damage fluences  J/cm^2

    if (strcmp(obj.material_stack,'SOI'))
        if (obj.verbose)
            display('    loading damage data for Si ...');
        end
        Fs = Fs*Si_factor/SiO2_factor;
        n = obj.nSi;
    else
        if (obj.verbose)        
            display('    loading damage data for Si3N4 ...');
        end
        Fs = Fs*Si3N4_factor/SiO2_factor;
        n = obj.nSi3N4;    
    end

    Es = sqrt(Fs/(1e-2)^2./10.^durations/1e-12*2/obj.c0/obj.e0/n);  % Peak field damage V/m

    % Return this inline function to give damage threshold (V/m) as fn of pulse
    % duration (ps)
    Ed = @(tau) interp1(durations,Es,log10(tau/1e-12));

end