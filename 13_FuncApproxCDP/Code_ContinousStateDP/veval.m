function v = veval(kprime, k, c, cheb, m)
   cons = k^m.alpha - kprime + (1-m.delta)*k;
   if (m.sigma==1)
       u = log(cons);
   else
       u = cons^(1-m.sigma)/(1-m.sigma);
   end
   v = -(u + m.beta*chebEval(kprime, c, cheb.order,cheb.lo, cheb.hi));
end
   
