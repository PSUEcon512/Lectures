%This is the function of the equilibrium constraints for the case
%with interstate borders in germany and a  a log distance cost
%
% PAUL GRIECO OCT.14.2013
function [c,ceq, dc, dceq] = model_constraints(rho,par,m)



% Read root
rho_ger = rho(1:m.num_pro_ger * (m.num_firms_ger + 1) ,1);    
rho_dnk = rho(m.num_pro_ger * (m.num_firms_ger + 1)+1 : m.num_pro_ger * (m.num_firms_ger + 1) + m.num_pro_dnk * (m.num_firms_dnk + 1) ,1) ;    
fe      = par(1 : m.num_firms_ger);         % firm fixed effect (note: all the firms are active in Germany)
betaLogD = par(m.num_firms_ger + 1);              % coefficient on distance 
betaFor = par(m.num_firms_ger + 2);
%betaState = par(m.num_firms_ger + 3);


% Germany
rho_ger_mat    = (reshape(rho_ger,m.num_firms_ger + 1,m.num_pro_ger))';
rho_ger_vec_f  = rho_ger_mat(:,1);
rho_ger_mat_nf = rho_ger_mat(:,2:m.num_firms_ger + 1);

%This line contains the cost function for the spline/intrastate estimation:
c_mat_nf_ger   = repmat(fe',m.num_pro_ger,1) ...
    + betaLogD .* m.dist_mat_ger ...
    + betaFor .* m.foreign_mat_ger;
%    + betaState.*m.state_mat_ger;


x1_ger         = log(rho_ger_mat_nf) - log(repmat(rho_ger_vec_f,1,m.num_firms_ger)) + c_mat_nf_ger + 1 ./(1-rho_ger_mat_nf);
x2_ger         = sum(rho_ger_mat_nf,2) + rho_ger_vec_f - 1;

temp_ger_mat   = [x1_ger , x2_ger];
temp_ger_vec   = reshape(temp_ger_mat',m.num_pro_ger * (m.num_firms_ger + 1),1);

%Denmark
rho_dnk_mat    = (reshape(rho_dnk,m.num_firms_dnk + 1,m.num_pro_dnk))';
rho_dnk_vec_f  = rho_dnk_mat(:,1);
rho_dnk_mat_nf = rho_dnk_mat(:,2:m.num_firms_dnk + 1);

%This line contains the cost function for the spline/intrastate estimation:
c_mat_nf_dnk   = repmat(fe(1:m.num_firms_dnk,1)',m.num_pro_dnk,1) ...
    + betaLogD .* m.dist_mat_dnk ...
    + betaFor .*  m.foreign_mat_dnk;
%    + betaState.* m.state_mat_dnk;  


x1_dnk         = log(rho_dnk_mat_nf) - log(repmat(rho_dnk_vec_f,1,m.num_firms_dnk)) + c_mat_nf_dnk + 1 ./(1-rho_dnk_mat_nf);
x2_dnk         = sum(rho_dnk_mat_nf,2) + rho_dnk_vec_f - 1;

temp_dnk_mat   = [x1_dnk , x2_dnk];
temp_dnk_vec   = reshape(temp_dnk_mat',m.num_pro_dnk * (m.num_firms_dnk + 1),1);


% Combine
ceq           = [temp_ger_vec ; temp_dnk_vec];
c             = [];

dc = [];
dceq = zeros(size(ceq,1));

%Step 1: Create Jacobian for German rhos....square
%matrix:(num_ger_firms+1)*num_pro_ger x same
%
%This is the first block of the full jacobian.  We create it by looping
%over the german firms, kroneckering the resulting calculations against a
%block template, and then summing the resulting "firm-level" matricies
%together.  
%
%r_0_matrix = sparse(kron(diag(rho_ger_vec_f), [ones(1,m.num_firms_ger),0;zeros(m.num_firms_ger,m.num_firms_ger+1)]));
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


%Step 3: Put it all together....Same as eqm_constraints but without
%parmeters...
dceq = [jac_ger zeros((m.num_firms_ger+1)*m.num_pro_ger, (m.num_firms_dnk+1)*m.num_pro_dnk); ...
        zeros((m.num_firms_dnk+1)*m.num_pro_dnk, (m.num_firms_ger+1)*m.num_pro_ger) jac_dnk; ];
