function z = priceFOC(p1)

global L c g;
p2 = p1';
%meshgrid outputs repeated rows, then repeated cols
[g2, g1] = meshgrid(g,g);

z = 1 - (((1+ exp(g2 - p2))./(1 + (exp(g1 - p1) + exp(g2 -p2)))).*(p1 - c));
