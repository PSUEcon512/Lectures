function out=solvelab(in, kcurrent, knext)

in=max(in,0.001);  % we do not want to program try inputs outside of (0,1) range
in=min(in,0.999);

global  mu delta alpha

C =  ( ( ( kcurrent^alpha)*(in^(1-alpha)) ) - knext + (1-delta)*kcurrent);

    if C>0; % if C is not positive we are in trouble, we set MUC very high
        MUC=1/C;
    else
        MUC=100000;
     end
              
     MUL=1/(1-in);             
     MPL = (1-alpha)*( kcurrent^alpha)*(in^(-alpha));
      
     out = (1-mu)*MUC*MPL - mu*MUL ;         