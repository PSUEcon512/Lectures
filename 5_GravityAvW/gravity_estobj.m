function F = gravity_estobj(beta,m)
temp     = beta(3)*(1-m.la) + beta(4)*(1-m.hi) + beta(5)*m.eu + beta(6)*m.na + beta(7)*m.as; 
ltc      = beta(2)*log(m.di) + temp;             % log of trade cost function  
F = log(m.tr) - log(m.eg) - log(m.ig) - beta(1) - (1-m.sigma)*ltc + (1-m.sigma)*log(m.Pex) + (1-m.sigma)*log(m.Pim) ;