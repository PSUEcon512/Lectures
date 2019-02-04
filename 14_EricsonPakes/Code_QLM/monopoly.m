%
% Quality ladder monopoly.
% Value and policy funciton iteration.
%

clear;
delete('monopoly.lst');
diary('monopoly.lst');

% Globals.
global L delta alpha c;

% Setup.
run setup;

% Product market equilibrium.
pi = zeros(L,1);
for i=1:L

    % Compute optimal prices. Solve FOCs.
    MaxIter = 400; TolFun = 1e-6; TolX = 1e-6; options = [1,TolFun,TolX,TolFun,MaxIter,1e-6];
    [pstar,info] = DogLeg('FOC',g(i),max(c,max(g)),options);
    if info(6)<=3
        disp('DogLeg terminated successfully.');
    else
        warning('There is a problem with DogLeg!');
        pause
    end

    % Demand.
    u = exp(g(i)-pstar);
    d = u./(1+u);

    % Profits.
    pi(i) = M.*d.*(pstar-c);

end

% Value function iteration.
disp(' ');
disp('Value function iteration');
disp('------------------------');
disp(' ');
V0 = zeros(L,1);
V1 = zeros(L,1);
x1 = zeros(L,1);
tic
iter = 0;
done = 0;
while ~done

    % Step 1: Maximization.
    for i=1:L

        % Compute optimal policy.
        if i==1
            x1(i) = (-1+sqrt(max(1,beta.*alpha.*(1-delta).*(V0(i+1)-V0(i)))))./alpha;
        elseif i==L
            x1(i) = (-1+sqrt(max(1,beta.*alpha.*delta.*(V0(i)-V0(i-1)))))./alpha;
        else
            x1(i) = (-1+sqrt(max(1,beta.*alpha.*((1-delta).*(V0(i+1)-V0(i))+delta.*(V0(i)-V0(i-1))))))./alpha;
        end

        % Update value function.
        p = transprob(i,x1(i));
        V1(i) = pi(i)-x1(i)+beta.*sum(V0.*p');

    end

    % Step 2: Convergence.
    tolV = max(abs(V1-V0));
    if tolV<tol
        done = 1;
    end
    V0 = V1;
    iter = iter+1;
    if (mod(iter,10)==0) | done
        disp(sprintf('iter=%d tolV=%g',iter,tolV));
    end

end
toc

% Print results.
disp(' ');
disp(' i         V(i)         x(i)');
disp('----------------------------');
for i=1:L
    disp(sprintf('%2d %12.8f %12.8f',i,V1(i),x1(i)));
end

% Policy function iteration.
disp(' ');
disp('Policy function iteration');
disp('-------------------------');
disp(' ');
V0 = zeros(L,1);
V1 = zeros(L,1);
x1 = zeros(L,1);
tic
iter = 0;
done = 0;
while ~done

    % Step 1: Maximization.
    P = zeros(L,1);
    Q = zeros(L,L);
    for i=1:L

        % Compute optimal policy.
        if i==1
            x1(i) = (-1+sqrt(max(1,beta.*alpha.*(1-delta).*(V0(i+1)-V0(i)))))./alpha;
        elseif i==L
            x1(i) = (-1+sqrt(max(1,beta.*alpha.*delta.*(V0(i)-V0(i-1)))))./alpha;
        else
            x1(i) = (-1+sqrt(max(1,beta.*alpha.*((1-delta).*(V0(i+1)-V0(i))+delta.*(V0(i)-V0(i-1))))))./alpha;
        end

        % Compute per-period profits and transition probabilities.
        P(i) = pi(i)-x1(i);
        Q(i,:) = transprob(i,x1(i));

    end

    % Step 2: Inversion.
    V1 = (eye(L,L)-beta.*Q)\P;

    % Step 3: Convergence.
    tolV = max(abs(V1-V0));
    if tolV<tol
        done = 1;
    end
    V0 = V1;
    iter = iter+1;
    if (mod(iter,1)==0) | done
        disp(sprintf('iter=%d tolV=%g',iter,tolV));
    end

end
toc

% Print results.
disp(' ');
disp(' i         V(i)         x(i)');
disp('----------------------------');
for i=1:L
    disp(sprintf('%2d %12.8f %12.8f',i,V1(i),x1(i)));
end

diary off;
