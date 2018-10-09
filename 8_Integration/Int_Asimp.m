function [ returnVal, N ] = Int_Asimp(f, a, b, tol )
%Int_Asimp: Adaptive quadrature with a tolerance, this could be much more
%efficient by re-using function evaluations. 

    assert(a < b);
    
    N = 10;
    oldV = Int_simp(f, a, b, N);
    maxN = 10000;
    
    while (N < maxN)
        N = 2*N;
        newV = Int_simp(f, a, b, N);
        %Check a relative tolerance:
        if (abs(newV - oldV)/(newV+1) < tol)
            break;
        end
        oldV = newV;
    end

    returnVal = newV; 
end
