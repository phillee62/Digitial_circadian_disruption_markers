function hr_phase_est = modify_hr_phase(hr_phase_est, t_phase)

    t_phase = t_phase - 1;
    
    if hr_phase_est > t_phase + 12
        hr_phase_est = hr_phase_est - 24;
    elseif hr_phase_est < t_phase - 12
        hr_phase_est = hr_phase_est + 24;
    end

end