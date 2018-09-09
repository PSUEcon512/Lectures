function pfval=gravity_priobj(P,m)

% m.tc = trade cost function, y=income 

pfval   = P.^(1-m.sigma) - ((m.tc).^(1-m.sigma) * ( m.y .* (P.^(m.sigma-1)) ));
%pfval(67) = (1-P(67))*100000000000;

