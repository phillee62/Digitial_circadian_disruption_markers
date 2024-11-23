function phase_std = find_phase_std(xM_phase)
    
    global X_std
        
    X = mvnrnd(xM_phase(:,1),xM_phase(:,2:4)*xM_phase(:,2:4).',10000);
    X_std = zeros(length(X),1); 
    mean_phase = measurement_function(xM_phase(:,1));
    
    for i=1:length(X) 
        X_std(i) = measurement_function(X(i,:)');   
        
        if X_std(i) > mean_phase + 12
            X_std(i) = X_std(i) - 24;
        elseif X_std(i) < mean_phase - 12
            X_std(i) = X_std(i) + 24;
        end
    end
    
    X_std = X_std(-12 <= X_std & X_std <= 36);
    phase_std = std(X_std);
    
end
