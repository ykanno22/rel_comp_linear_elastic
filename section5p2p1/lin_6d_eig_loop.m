%
clear;
close all;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Flag.save = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
if Flag.save == 1
    delete('lin_6d_eig_loop.log');
    diary( 'lin_6d_eig_loop.log');
end
%
param.eps  = 0.10;
param.delta = 0.10;
%
num_Iter = 10;
Length = 60.0; % Length of the model in [mm]
Width  = 20.0; % Width               in [mm]
NXE = 30;       % Number of rows in the x direction
NYE = 10;       % Number of rows in the y direction
%
load('lin_isotropic_data_set.mat');
num.data = size(list_eps, 2);
num.dim   = 6;
num.dim_x = num.dim / 2;
num.variable = num.data * (num.dim + 1);
%
data_x(1:3,:) = list_eps;
data_x(4:6,:) = list_sigma;
clear list_eps list_sigma
% <---
% Some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare matrices & vectors
% --->
vec_h = cell(1,num.data);
for j=1:num.data
    vec_h{j} = [data_x(:,j); -1];
end
%
matH  = zeros(num.data, num.dim+1);
for j=1:num.data
    matH(j,:) = vec_h{j}';
end
matQ = matH' * matH;
% <---
% Prepare matrices & vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning via eigenvalue analysis
% --->
[eig_vec, eig_val] = eig(matQ);
eig_val = diag(eig_val);
[eig_val, idx_sort] = sort(eig_val);
eig_vec = eig_vec(:,idx_sort);
%
fprintf(' ============================================= \n');
fprintf('   eigenvalue     square_error \n');
for i=1:7
    fprintf('     %3.5d    %3.5d \n',...
        eig_val(i),...
        eig_vec(:,i)' * (matH') * matH * eig_vec(:,i));
end
fprintf(' ============================================= \n');
%
vec_sol = cell(1, num.dim);
for k=1:3
    vec_sol{k} = eig_vec(:,k);
end
%
vec_a  = zeros(6,3);
cons_c = zeros(1,3);
for k=1:3
    vec_a(:,k)  = vec_sol{k}(1:6);
    cons_c(:,k) = vec_sol{k}(7);
end
% <---
% Learning via eigenvalue analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num.plane = 3;
np = num.plane;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #points included in uncertainty set
% --->
num_confi = num.data;
cur_delta = (1-param.eps)^num.data;

if cur_delta > param.delta
    fprintf(' confidence level "delta" is too small \n');
    return;
end

while cur_delta <= param.delta
    num_confi = num_confi - 1;
    cur_delta = cur_delta + ...
        ( nchoosek(num.data,num_confi)...
        * ( (1-param.eps)^num_confi )...
        * ( param.eps^(num.data-num_confi) ) );
end

fprintf(' ============================================= \n');
fprintf('   "1 - %4.3f reliability" with "1 - %4.3f confidence" \n',...
    param.eps, param.delta );
fprintf('   (#conf-1)/(#data) = %g/%g: LHS = %4.3e >  delta = %4.3e \n',...
    num_confi, num.data, cur_delta, param.delta );
num_confi = num_confi + 1;
cur_delta = 0;
for jj=num_confi:num.data
    cur_delta = cur_delta +...
        ( nchoosek(num.data,jj)...
        * ( (1-param.eps)^jj )...
        * ( param.eps^(num.data-jj) ) );
end
fprintf('   (#conf)  /(#data) = %g/%g: LHS = %4.3e <= delta = %4.3e \n',...
    num_confi, num.data, cur_delta, param.delta );
fprintf(' ============================================= \n');
% <---
% #points included in uncertainty set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Order statistics
% --->
order_s =...
    max(abs( (vec_a' * data_x) - cons_c' ));
[~,Idx_sort] = sort(order_s);
val_tau = order_s(Idx_sort(num_confi));
% <---
% Order statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To change the size of the mesh, alter the next statements
% %%%%% Introduction to Finite Element Analysis Using MATLAB and Abaqus
% --->
thick = 5.0;      % Tickness in mm
E  = 20.0;      % Elastic modulus in  100*MPa
nu = 0.3;      % Poisson's ratio 
E  = mean_Young;
nu = mean_Poisson;
%
% Force = 200.0;    % external force in kN
Force = 0.25 / NXE; %%% distributing: Force / (Length * thhick) in [kN/mm2]
%
dhx = Length/NXE; % Element size in the x direction
dhy = Width/NYE;  % Element size in the x direction
X_origin = 0.0;       % X origin of the global coordinate system
Y_origin = 0.0;   % Y origin of the global coordinate system
% <--
% To change the size of the mesh, alter the next statements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the mesh 
% --->
nne = 4;
nodof = 2;
eldof = nne*nodof;
%
[nnd, nel, geom, connec] = Q4_mesh_fun(NXE, NYE, X_origin, Y_origin, dhx, dhy); 
% <--
% Generate the mesh 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form the elastic matrix for plane stress 
% --->
dee = formdsig(E, nu);
% <--
% Form the elastic matrix for plane stress 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions
% --->
%%%% Restrain in all directions the nodes situated @ (x = 0)
nf = ones(nnd, nodof);
%
fixed_dig = (geom(:,1)==0);
nf(fixed_dig,:) = 0;
clear fixed_dig
%
nd = sum(sum(nf));
idx_bottom_right_ydirec = nd - (2 * NYE);
%
ii = 0;
for i=1:nnd
    for j=1:nodof
        if nf(i,j) ~= 0
            ii = ii + 1;
            nf(i,j) = ii;
        end
    end
end

%%%% Apply uniform distributed load at the nodes at y=Width
Nodal_loads= zeros(nnd, 2);
Ind_top_nodes = find(geom(:,2) == Width);
Nodal_loads(Ind_top_nodes,2)      = -Force;
Nodal_loads(Ind_top_nodes(end),2) = -Force/2;
%
Nodal_loads = sparse(Nodal_loads);
% <--
% Boundary conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble the global force vector
% --->
vec_f = zeros(nd,1);
%
Idx_tmp = find(nf(:,1));
vec_f(nf(Idx_tmp,1)) = Nodal_loads(Idx_tmp,1);
clear Idx_tmp
%
Idx_tmp = find(nf(:,2));
vec_f(nf(Idx_tmp,2)) = Nodal_loads(Idx_tmp,2);
clear Idx_tmp
%
vec_f = sparse(vec_f) * param_scale;
% <---
% Assemble the global force vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Form the matrix containing the abscissas and the weights of Gauss points
% --->
ngp  = 2;
samp = gauss(ngp);
% <---
% Form the matrix containing the abscissas and the weights of Gauss points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical integration and assembly of the global stiffness matrix
% --->
mat_K = sparse(nd, nd);
mat_B = [];
integ_weight = [];
for i=1:nel
    %%%% coordinates of the nodes of element i, & its steering vector -->
    [coord, g] = elem_q4_fun(i, nne, nodof, geom, connec, nf);
    ke = sparse(eldof, eldof);
    for ig=1:ngp
        wi = samp(ig,2);
        for jg=1:ngp
            wj = samp(jg, 2);
            %%%% Derivative of shape functions in local coordinates -->
            [der, fun] = fmlin(samp, ig, jg);
            %%%% Compute Jacobian matrix -->
            jac = der * coord;
            %%%% Compute determinant of Jacobian matrix -->
            d = det(jac);
            %%%% Derivative of shape functions in global coordinates -->
            deriv = jac \ der;
            %%%% Form matrix [B] -->
            bee = formbee(deriv, nne, eldof);
            %%%% Integrate stiffness matrix -->
            integ_weight = [integ_weight, d * thick * wi * wj];
            mat_Be = zeros(3, nd);
            for jj=find(g)
                mat_Be(:,g(jj)) = bee(:,jj);
            end
            mat_B = [mat_B; mat_Be];
        end
    end
end
mat_K = mat_B' * kron(diag(integ_weight), dee) * mat_B;
% <---
% Numerical integration and assembly of the global stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conventional FEM
% --->
model_u = mat_K \ vec_f;

fprintf(' ============================================= \n');
fprintf(' size of elastic body = %3.2f mm x %3.2f mm \n',...
    Length, Width);
fprintf('     mesh size = %g x %g \n',...
    NXE, NYE);
fprintf(' conventional FEM \n');
fprintf('     Young = %3.5f MPa ;  Poisson = %3.5f \n',...
    mean_Young * param_scale, mean_Poisson);
fprintf('     displacement = %3.5f mm \n',...
    model_u(idx_bottom_right_ydirec) );
fprintf('          (right bottom node, y-direction) \n');
fprintf(' data-driven solver \n');
fprintf('     #data-points = %g \n', num.data);
% <---
% Conventional FEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% variables : [eps]    \in (3*ns)                   %%%%
%%%% variables : [sig]    \in (3*ns)                   %%%%
%%%% variables : [u]      \in (nd)                     %%%%
nm = nel;
ns = (ngp*ngp) * nm;  num.stress = ns;  num.strain = ns;

mat_N = (kron(diag(integ_weight), eye(3)) * mat_B)';
mat_cons_c = repmat(cons_c', 1, ns);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LP for computing bounds
% --->
his_load_fact = zeros(1, num_Iter+1);
his_model_u = zeros(nd, num_Iter+1);
his_min_u  = zeros(nd, num_Iter+1);
his_max_u  = zeros(nd, num_Iter+1);

incr_load_factor = 1 / num_Iter;

for Iter = 0:num_Iter
    fprintf('     ----------------------------------------- \n');
    load_fact = Iter * incr_load_factor;
    
    vec_f_cur = load_fact * vec_f;
    model_u = mat_K \ vec_f_cur;
    
    cvx_begin quiet
        cvx_solver sedumi
        cvx_precision best
        variable min_u(nd,1)
        variable min_eps(3,ns)
        variable min_sig(3,ns)
        minimize( min_u(idx_bottom_right_ydirec) )
        subject to
            reshape(min_eps, 3*ns, 1) == mat_B * min_u;
            mat_N * reshape(min_sig, 3*ns, 1) == vec_f_cur;
            max(max( abs( (vec_a' * [min_eps; min_sig]) - mat_cons_c ))) <= val_tau;
    cvx_end
    cvx_status_min = cvx_status;
    
    cvx_begin quiet
        cvx_solver sedumi
        cvx_precision best
        variable max_u(nd,1)
        variable max_eps(3,ns)
        variable max_sig(3,ns)
        maximize( max_u(idx_bottom_right_ydirec) )
        subject to
            reshape(max_eps, 3*ns, 1) == mat_B * max_u;
            mat_N * reshape(max_sig, 3*ns, 1) == vec_f_cur;
            max(max( abs( (vec_a' * [max_eps; max_sig]) - mat_cons_c ))) <= val_tau;
    cvx_end
    
    fprintf('     load_factor = %2.1f : FEM disp. = %3.5f mm \n',...
        load_fact, model_u(idx_bottom_right_ydirec) );
    fprintf('       displacement in [%3.5f, %3.5f] mm \n',...
        min_u(idx_bottom_right_ydirec), max_u(idx_bottom_right_ydirec) );
    fprintf('       cvx_status = %s, %s \n', cvx_status_min, cvx_status );
    
    his_model_u(:,Iter+1)  = model_u;
    his_min_u(:,Iter+1)  = min_u;
    his_max_u(:,Iter+1)  = max_u;
    his_load_fact(1,Iter+1) = load_fact;
end
% <---
% LP for computing bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' ============================================= \n');
% <---
% Output the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw figures
% --->
figure;
plot(his_model_u(idx_bottom_right_ydirec,:), his_load_fact, 'bs:',...
    'LineWidth', 1.0, 'MarkerSize',6);
hold on;
plot(his_max_u(idx_bottom_right_ydirec,:), his_load_fact, 'rv-',...
    'LineWidth', 1.0, 'MarkerSize',6);
plot(his_min_u(idx_bottom_right_ydirec,:), his_load_fact, 'r^-',...
    'LineWidth', 1.0, 'MarkerSize',6);
xlabel('Displacement (mm)', 'Interpreter', 'latex');
ylabel('Load factor, $\lambda$', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
xlim([1.1*min(his_model_u(idx_bottom_right_ydirec,:)), 1.0]);
ylim([0, 1.1]);
if Flag.save == 1
     print('-depsc2', '-painters', 'lin_6d_eig_loop');
end
% <---
% Draw figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Flag.save == 1
    num
    save('lin_6d_eig_loop.mat');
end

diary off;

