function [f grad_f] = likelihood_func(X,m)

rho_ger = X(1: m.num_pro_ger*(m.num_firms_ger + 1),1);
rho_dnk = X(m.num_pro_ger*(m.num_firms_ger + 1)+1 : m.num_pro_ger*(m.num_firms_ger + 1) + m.num_pro_dnk*(m.num_firms_dnk + 1),1);

part1 = log((reshape(rho_ger,m.num_firms_ger + 1,m.num_pro_ger))') .* m.y_matrix_ger ;  % num_pro_ger  x  (num_firms_ger + 1) matrix
part2 = log((reshape(rho_dnk,m.num_firms_dnk + 1,m.num_pro_dnk))') .* m.y_matrix_dnk ;

f  = -(sum(sum(part1)) + sum(sum(part2)));   % "-" since knitro minimizes objective function 



grad_f = zeros(size(X,1),1);
grad_f(1:(m.num_firms_ger + 1)*m.num_pro_ger,1) = - ( 1./rho_ger ) .* reshape(m.y_matrix_ger',(m.num_firms_ger + 1)*m.num_pro_ger,1) ;
grad_f((m.num_firms_ger + 1)*m.num_pro_ger+1:(m.num_firms_ger + 1)*m.num_pro_ger+(m.num_firms_dnk + 1)*m.num_pro_dnk,1) = - ( 1./rho_dnk ) .* reshape(m.y_matrix_dnk',(m.num_firms_dnk + 1)*m.num_pro_dnk,1) ;


