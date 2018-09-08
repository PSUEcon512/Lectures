function [fval fjac] = cournot(q)
    % Set parameters:
    c = [0.6; 0.8];
    eta = 1.6;
    e = -1/eta; %For readability

    %Evaluate function: 
    fval = sum(q)^e + e*sum(q)^(e-1)*q - diag(c)*q;

    %Evaluate Jacobian: 
    fjac = e*sum(q)^(e-1)*ones(2,2) + ...
           e*sum(q)^(e-1)*eye(2) + ...
           (e-1)*e*sum(q)^(e-2)*q*[1 1] - diag(c);
end