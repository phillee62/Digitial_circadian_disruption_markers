function dxMdt = dxMdt(tspan, y, xM, k_noise, light)  
    x_and_M = reshape(y, size(xM));
    x_mean = x_and_M(:,1);
    dim = numel(x_mean);
    M = x_and_M(:,2:(dim+1));
    dx_meandt = zeros(dim,1); 
    dMdt = zeros(dim, dim);  
    
    %% Find the average velocity (velocity of the Gaussian ellipsoid mean)
    for i = 1:dim
        x_plus_span1 = x_mean + M(:,i);
        x_plus_span2 = x_mean - M(:,i);
        v_mean_plus_span1 = simpler_circadian(tspan, x_plus_span1, light);
        v_mean_plus_span2 = simpler_circadian(tspan, x_plus_span2, light);
        dx_meandt = dx_meandt + 0.5.*(1./dim).*(v_mean_plus_span1 + v_mean_plus_span2);  
    end
    
    %% Find the velocity of level set
    for i = 1:dim
        e = zeros(dim,1);
        x_plus_span = x_mean + M(:,i);
        xi_velocity = simpler_circadian(tspan, x_plus_span, light);
        e(i) = 1;
        dxidt = xi_velocity - dx_meandt + 0.5.*k_noise.*inv(M.')*e;
        dMdt(:,i) = dxidt;
    end     
    
    %% Return ODE output
    dxMdt = [dx_meandt, dMdt]; 
    dxMdt = reshape(dxMdt, [numel(dxMdt),1]);
    
end

