function [] = plot_pulse_helper(obj)
    [Nprobe,~] = size(obj.A_t_list);
    figure; clf; hold all;
    for i = (1:Nprobe)
        plot(obj.ts./1e12,real(obj.A_t_list(i,:))./1e6);
    end
    xlabel('time (ps)');
    ylabel('real (E(t)) (MV/m)');

    figure; clf; hold all;            
    for i = (1:Nprobe)
        plot(obj.fs./1e12,1e-9*abs(fftshift(fft(obj.A_t_list(i,:)))));
    end
    xlabel('(frequency - f0) (THz)');
    ylabel('abs (E(f)) (GV/m/Hz)');
    xlim([-50 50])
end