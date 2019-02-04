function W = getW(Vin, x)
global L M c alpha delta beta g CRIT;

%Construct Continuation Value as a three-dimensional array:
%(Player 1 State, Player 2 State, Player 1 Shift)
%The third dimension conditions on only the ``shift'' in player 1's state
%between t and t+1, i.e., does his quality improve, stay the same, or drop.
%

x2 = x';
%x2 = x;
cont = zeros(L,L);
V = [zeros(1,L); Vin; zeros(1,L)];
for i=-1:1
    shft = i+2;
	rws = 2+i:L+1+i;
    
    cont(:, 1, shft) = V(rws, 1).*((1+delta*alpha.*x2(:,1))./(1 + alpha.*x2(:,1))) + ...
        V(rws,2).*(((1-delta)*alpha.*x2(:,1))./(1 + alpha.*x2(:,1)));
   
    cont(:,L,shft) = V(rws, L).*((1-delta+alpha.*x2(:,L))./(1 + alpha.*x2(:,L))) + ...
        V(rws, L-1).*(delta./(1 + alpha*x2(:,L)));

    cont(:, 2:L-1,shft) = V(rws, 3:L).*((1-delta)*alpha.*x2(:, 2:L-1))./(1 + alpha.*x2(:, 2:L-1)) + ...
        V(rws, 2:L-1).*(1-delta+delta*alpha.*x2(:, 2:L-1))./(1 + alpha.*x2(:, 2:L-1)) + ...
        V(rws, 1:L-2).*(delta./(1 + alpha*x2(:, 2:L-1)));
end
W = cont;