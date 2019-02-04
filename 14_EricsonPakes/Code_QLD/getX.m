function x = getX(W);
global L M c alpha delta beta g CRIT;

%Construct Player 1's optimal strategy in each state. Recall that this
%depends only on the difference in continuation values for player 1's next
%period state, which is why we first construct fdWb (first-difference W
%backward) and fdWf (first difference W backward). Which reflect the effect
%of improvements and falls in player 1's state, conditioning on the current
%state.

%fdW = W(2:L, :) - W(1:L-1, :);
fdWb = W(:,:,2) - W(:,:,1); 
fdWf = W(:,:,3) - W(:,:,2);
x(1,:) = (-1 + sqrt(beta*alpha*(1-delta).*fdWf(1,:)))./alpha;
x(L,:) = (-1 + sqrt(beta*alpha*delta*fdWb(L,:)))./alpha;
x(2:L-1,:) = (-1 + sqrt(beta*alpha*((1-delta)*fdWf(2:L-1,:)+ delta*fdWb(2:L-1,:))))./alpha;
x = max(zeros(L,L), x.*(imag(x)==0));