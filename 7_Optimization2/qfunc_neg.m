function [ f, grad ] = qfunc_neg( x )
% A quartic function to minimize, complete with gradient:
% Recall we want to minimize:
%
% $$f(x)=(x_1-x_2)^4+2x_1^2+x_2^2-x_1+2x_2$$
%
  f = -(x(1) - x(2))^4 + 2*x(1)^2 + x(2)^2 - x(1) + 2*x(2);
% The Gradient is,  
%
% $$\nabla f(x) = \left( \begin{array}{c} 4(x_1-x_2)^3+4x_1-1 \\ -4(x_1-x_2)^3+2x_2+2  \end{array} \right) $$
%
  grad = -[4*(x(1)-x(2))^3+4*x(1)-1;...
          -4*(x(1)-x(2))^3+2*x(2)+2];  
      
%Won't need the Hessian since we are using BFGS...
end

