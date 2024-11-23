function dydt = simpler_circadian(t,y,light)
    x = y(1); xc = y(2); n = y(3);
    
    %% Define parameters
    try 
        I = light(mod(floor(t*60/5),length(light))+1);
    catch
        I = light(end);
    end
    
    alp0 = 0.05; bet = 0.0075; G = 33.75; p = 0.6; k = 0.55; I0 = 9500; tau = 24.2; mu = 0.23;

    alpha = alp(I, I0, alp0, p); 
    Bhat = G.*(1-n).*alpha;
    B = Bhat.*(1-0.4.*x).*(1-0.4.*xc);
    
    %% ODEs
    dxdt = (pi/12).*(xc+B);
    dxcdt = (pi/12).*(mu.*(xc-(4.*(xc.^3))./3)-x.*((24./(0.99669*tau)).^2+k.*B));
    dndt = 60.*(alpha.*(1-n)-bet.*n);
    dydt = [dxdt; dxcdt; dndt];
    
    function alpha = alp(I, I0, alp0, p)
        alpha = alp0.*(I./I0).^p;
    end
end

