%This is the function of the equilibrium constraints for the case
%with interstate borders in germany and a log distance cost
%
% PAUL GRIECO OCT.14.2013
function [c,ceq, dc, dceq] = eqm_constraints(X,m)




% Read root
rho_ger = X(1:m.num_pro_ger * (m.num_firms_ger + 1) ,1);    
rho_dnk = X(m.num_pro_ger * (m.num_firms_ger + 1)+1 : m.num_pro_ger * (m.num_firms_ger + 1) + m.num_pro_dnk * (m.num_firms_dnk + 1) ,1) ;    
fe      = X(m.num_pro_ger * (m.num_firms_ger + 1) + m.num_pro_dnk * (m.num_firms_dnk + 1) +1 : m.num_pro_ger * (m.num_firms_ger + 1) + m.num_pro_dnk * (m.num_firms_dnk + 1) + m.num_firms_ger);         % firm fixed effect (note: all the firms are active in Germany)
betaLogD   = X(m.num_pro_ger * (m.num_firms_ger + 1) + m.num_pro_dnk * (m.num_firms_dnk + 1) + m.num_firms_ger + 1, 1);              % coefficient on distance 
% coefficient of the International and Interstate borders
betaFor   = X(m.num_pro_ger * (m.num_firms_ger + 1) + m.num_pro_dnk * (m.num_firms_dnk + 1) + m.num_firms_ger + 2, 1); 
%betaState = X(m.num_pro_ger * (m.num_firms_ger + 1) + m.num_pro_dnk * (m.num_firms_dnk + 1) + m.num_firms_ger + 3, 1); 


% Germany
rho_ger_mat    = (reshape(rho_ger,m.num_firms_ger + 1,m.num_pro_ger))';
rho_ger_vec_f  = rho_ger_mat(:,1);
rho_ger_mat_nf = rho_ger_mat(:,2:m.num_firms_ger + 1);

c_mat_nf_ger   = repmat(fe',m.num_pro_ger,1) ...
    + betaLogD .* m.dist_mat_ger ...
    + betaFor .* m.foreign_mat_ger;
%    + betaState.*m.state_mat_ger;
%Here I add the interstate border effect and the radius effect
                                                                 

x1_ger         = log(rho_ger_mat_nf) - log(repmat(rho_ger_vec_f,1,m.num_firms_ger)) + c_mat_nf_ger + 1 ./(1-rho_ger_mat_nf);
x2_ger         = sum(rho_ger_mat_nf,2) + rho_ger_vec_f - 1;

temp_ger_mat   = [x1_ger , x2_ger];
temp_ger_vec   = reshape(temp_ger_mat',m.num_pro_ger * (m.num_firms_ger + 1),1);

%Denmark the code for the danish market need not be changed for the
%presence of german interstate border effect
rho_dnk_mat    = (reshape(rho_dnk,m.num_firms_dnk + 1,m.num_pro_dnk))';
rho_dnk_vec_f  = rho_dnk_mat(:,1);
rho_dnk_mat_nf = rho_dnk_mat(:,2:m.num_firms_dnk + 1);

c_mat_nf_dnk   = repmat(fe(1:m.num_firms_dnk,1)',m.num_pro_dnk,1) ...
    + betaLogD .* m.dist_mat_dnk ...
    + betaFor .* m.foreign_mat_dnk;
%    + betaState.*m.state_mat_dnk;  

x1_dnk         = log(rho_dnk_mat_nf) - log(repmat(rho_dnk_vec_f,1,m.num_firms_dnk)) + c_mat_nf_dnk + 1 ./(1-rho_dnk_mat_nf);
x2_dnk         = sum(rho_dnk_mat_nf,2) + rho_dnk_vec_f - 1;

temp_dnk_mat   = [x1_dnk , x2_dnk];
temp_dnk_vec   = reshape(temp_dnk_mat',m.num_pro_dnk * (m.num_firms_dnk + 1),1);


% Combine
ceq           = [temp_ger_vec ; temp_dnk_vec];
c             = [];

dc = [];
dceq = zeros(size(X,1), size(ceq,1));

%Step 1: Create Jacobian for German rhos....square
%matrix:(num_ger_firms+1)*num_pro_ger x same
%
%This is the first block of the full jacobian.  We create it by looping
%over the german firms, kroneckering the resulting calculations against a
%block template, and then summing the resulting "firm-level" matricies
%together.  
%
%r_0_matrix = sparse(kron(diag(rho_ger_vec_f), [ones(1,m.num_firms_ger),0;zeros(m.num_firms_ger,m.num_firms_ger+1)]));
%This code does not have to be altered for the interstate border case
r_0_matrix = kron(-diag(1./rho_ger_vec_f), [ones(1,m.num_firms_ger),0;zeros(m.num_firms_ger,m.num_firms_ger+1)]);

jac_ger = r_0_matrix;

for j =1:m.num_firms_ger
    temp = zeros(m.num_firms_ger+1,m.num_firms_ger+1);
    temp(1+j,j) = 1;
 %   jac_j = sparse(kron(diag(rho_ger_mat_nf(:,j)),temp));
    %Next line computes the derivatives, rest is just getting them into the
    %correct location in the matrix.
    DcDr = (1./rho_ger_mat_nf(:,j)) + (1-rho_ger_mat_nf(:,j)).^(-2);
    jac_j = kron(diag(DcDr),temp);
    jac_ger = jac_ger + jac_j;
end

add_up_ger = [zeros(1,m.num_firms_ger) 1; 
              zeros(m.num_firms_ger, m.num_firms_ger) ones(m.num_firms_ger,1)];
%full_aug = sparse(kron(eye(m.num_pro_ger) ,add_up_ger));
full_aug = kron(eye(m.num_pro_ger) ,add_up_ger);
jac_ger = jac_ger + full_aug;




%Step 2: Create Jacobian for Danish market rhos..again a square matrix:
% (num_firms_dnk+1)*num_pro_dnk x same
r_0_matrix = kron(-diag(1./rho_dnk_vec_f), [ones(1,m.num_firms_dnk),0;zeros(m.num_firms_dnk,m.num_firms_dnk+1)]);

jac_dnk = r_0_matrix;

for j =1:m.num_firms_dnk
    temp = zeros(m.num_firms_dnk+1,m.num_firms_dnk+1);
    temp(1+j,j) = 1;
 %   jac_j = sparse(kron(diag(rho_ger_mat_nf(:,j)),temp));
    %Next line computes the derivatives, rest is just getting them into the
    %correct location in the matrix.
    DcDr = (1./rho_dnk_mat_nf(:,j)) + (1-rho_dnk_mat_nf(:,j)).^(-2);
    jac_j = kron(diag(DcDr),temp);
    jac_dnk = jac_dnk + jac_j;
end

add_up_dnk = [zeros(1,m.num_firms_dnk) 1; 
              zeros(m.num_firms_dnk, m.num_firms_dnk) ones(m.num_firms_dnk,1)];
%full_aug = sparse(kron(eye(m.num_pro_ger) ,add_up_ger));
full_aug = kron(eye(m.num_pro_dnk) ,add_up_dnk);
jac_dnk = jac_dnk + full_aug;


%Step 3: Fixed Effects
%These will cover all columns for the fixed effects rows.
%Need not be changed for the interstate border effects
FE_GER = [eye(m.num_firms_ger) zeros(m.num_firms_ger,1)];
FE_DNK = [ [eye(m.num_firms_dnk) zeros(m.num_firms_dnk, 1) ]; zeros(m.num_firms_ger - m.num_firms_dnk, m.num_firms_dnk+1)];
%Rows of Jacobian pertaining to fixed effects...
J_FE = [repmat(FE_GER, 1, m.num_pro_ger) repmat(FE_DNK, 1, m.num_pro_dnk)];


%Step 4: Distance Parameter - Affects all of the non-fringe (or non-adding
%up) constraints. Derivative is just the distance to a project -1...
%Notice we transpose in both lines to get a row vector taken row-wise
BETALD_GER = ([m.dist_mat_ger zeros(m.num_pro_ger,1)])';
BETALD_GER = BETALD_GER(:)';
BETALD_DNK = ([ m.dist_mat_dnk  zeros(m.num_pro_dnk,1)])';
BETALD_DNK = BETALD_DNK(:)';
J_BETALD = [BETALD_GER BETALD_DNK];

%Step 5: International Border Dummy
BETAFOR_GER = [ones(1,m.num_for_ger) zeros(1,m.num_firms_ger-m.num_for_ger+1)];
BETAFOR_DNK = [zeros(1,m.num_firms_dnk-m.num_for_dnk) ones(1,m.num_for_dnk) 0];
J_BETAFOR = [repmat(BETAFOR_GER, 1, m.num_pro_ger) repmat(BETAFOR_DNK, 1, m.num_pro_dnk)];

%Step 5.1: Interstate german dummies
%BETASTATE_GER = ([m.state_mat_ger zeros(m.num_pro_ger,1)])';
%BETASTATE_GER = BETASTATE_GER(:)';
%BETASTATE_DNK = ([m.state_mat_dnk zeros(m.num_pro_dnk,1)])';
%BETASTATE_DNK = BETASTATE_DNK(:)';
%J_BETASTATE = [BETASTATE_GER BETASTATE_DNK];


%Step 6: Put it all together....
dceq = [jac_ger zeros((m.num_firms_ger+1)*m.num_pro_ger, (m.num_firms_dnk+1)*m.num_pro_dnk); ...
        zeros((m.num_firms_dnk+1)*m.num_pro_dnk, (m.num_firms_ger+1)*m.num_pro_ger) jac_dnk; ...
        J_FE; J_BETALD; J_BETAFOR];







