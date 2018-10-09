function [ ret] = FindNormPointsWeights(wx)
% System of Equations to solve for nodes and weights of 3 point Gaussian
% Quadrature: 
  % Moments 0-5 of the normal distribution are:
  mom = [1, 0, 1, 0, 3, 0]';
  
  %w is the weights, x are the nodes
  w = wx(1:3); 
  x = wx(4:end)';
  
  % Discrete approximation of moments 0-5 using these nodes and weights.
  rhs = [ w*ones(3,1), w*x, w*(x.^2), w*(x.^3), w*(x.^4), w*(x.^5) ]';
  
  %Equations we want to solve:
  ret = mom - rhs;
end
