function [x , fx] = bisection(f, a, b, tolx, tolf) 
  % Make sure conditions for bisection apply,
  assert(b > a);
  assert(f(a) < 0);
  assert(f(b) > 0); 
  
  % MATLAB's way of allowing default parameters: 
  assert(nargin >= 3)
  if (nargin == 3)
      tolx = 1e-15; %This tolerance is to let program stop when bracket is small
      tolf = 1e-15; %This tolerance is to stop if it happens to compute a root
  end
  
  % Finally, the actual algorithm:
  fc = 1;
  while ((b - a) > tolx) && (abs(fc) > tolf)
      c = (a+b)/2;
      fc = f(c);
      if (fc < 0)
          a = c;
      else
          b = c;
      end
  end
  
  %Set return values
  x = c;
  fx = fc;
end
  