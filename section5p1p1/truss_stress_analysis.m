%
clear
close all
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Flag.save = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
if Flag.save == 1
    delete('prog2_truss_stress_analysis.log');
    diary( 'prog2_truss_stress_analysis.log');
end
%
param.eps  = 0.10;
param.delta = 0.10;
%
load('data_set.mat');
%
num.data = length(samples_eps);
num.dim  = 2;
data_x(1,:) = samples_eps';
data_x(2,:) = samples_sig';
clear samples_eps samples_sig
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
% Manifold-extracting via eigenvalue analysis
% --->
[eig_vec, eig_val] = eig(matQ);
eig_val = diag(eig_val);
[eig_val, idx_sort] = sort(eig_val);
eig_vec = eig_vec(:,idx_sort);
%
par_a = eig_vec(1, 1);
par_b = eig_vec(2, 1);
par_c = eig_vec(3, 1);
fprintf(' ============================================= \n');
fprintf('   eigenvalues = \n');
fprintf('        %3.5d \n',...
    eig_val);
% <---
% Manifold-extracting via eigenvalue analysis
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
    abs( (par_a * data_x(1,:))...
    + (par_b * data_x(2,:)) - par_c);
[~,Idx_sort] = sort(order_s);
val_tau = order_s(Idx_sort(num_confi));
% <---
% Order statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare figures: manifolds
% --->
x_plot = (min(data_x(1,:)) : max(data_x(1,:))/10 :max(data_x(1,:)));
x_plot = [x_plot, max(data_x(1,:))];
nx_plot = size(x_plot,2);
y_plot = zeros(3,nx_plot);
for i=1:nx_plot
    ss = x_plot(i);
    y_plot(1,i) = (par_c - (par_a * ss)) / par_b;
    y_plot(2,i) = ((par_c + val_tau) - (par_a * ss)) / par_b;
    y_plot(3,i) = ((par_c - val_tau) - (par_a * ss)) / par_b;
    clear ss
end
clear nx_plot
% <---
% Prepare figures: manifolds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truss data
% --->
[dll,matBt,coord_x,ir,irr,ird] = member(2,1);
%
nk = size(coord_x,1);  num.node   = nk;
nd = size(matBt,1);    num.degree = nd;
nm = size(matBt,2);    num.member = nm;
%
vec_cs = 10.0 * ones(nm,1);
%
% [~] = draw_cs(coord_x, irr, vec_cs);
% if Flag.save == 1
%     print('-depsc2', '-painters', '10bar_truss');
% end
% <---
% Truss data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vector
% --->
vec_p = zeros(nd,1);
vec_p(2) = -15.0;
vec_p(4) = -15.0;
% <---
% Load vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffness matrix
% --->
matK = matBt * sparse(diag(mean_Young * vec_cs ./ dll)) * matBt';
ref_u = matK \ vec_p;
matL  = diag(1 ./ dll) * matBt';
matN  = matBt * diag(vec_cs);
ref_eps = matL * ref_u;
ref_sig = mean_Young * ref_eps;
% <---
% Stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LP for computing bounds
% --->
min_sig_all = zeros(nm,1);
max_sig_all = zeros(nm,1);
for iMember = 1:nm
    cvx_begin quiet
        cvx_solver sedumi
        cvx_precision best
        variable min_u(nd,1)
        variable min_eps(nm,1)
        variable min_sig(nm,1)
        minimize( min_sig(iMember) )
        subject to
            min_eps == matL * min_u;
            matN * min_sig == vec_p;
            abs( (par_a * min_eps) + (par_b * min_sig) - par_c ) <= val_tau;
    cvx_end
    fprintf('       cvx_status = %s \n', cvx_status );

    cvx_begin quiet
        cvx_solver sedumi
        cvx_precision best
        variable max_u(nd,1)
        variable max_eps(nm,1)
        variable max_sig(nm,1)
        maximize( max_sig(iMember) )
        subject to
            max_eps == matL * max_u;
            matN * max_sig == vec_p;
            abs( (par_a * max_eps) + (par_b * max_sig) - par_c ) <= val_tau;
    cvx_end
    fprintf('       cvx_status = %s \n', cvx_status );
    
    min_sig_all(iMember) = min_sig(iMember);
    max_sig_all(iMember) = max_sig(iMember);
end
% <---
% LP for computing bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('   #members = %g ; cs(1) = %2.2f [cm^2] \n',...
    nm, nd);
fprintf('   dll(1) = %2.2f [m] ; cs(1) = %2.2f [cm^2] \n',...
    dll(1), vec_cs(1));
fprintf('   load = %2.2f [kN] \n',...
    norm(vec_p,Inf) / 10);
for iMember=1:nm
    fprintf('   ref_sig = %2.4e in [%2.4e, %2.4e]  * 10^6 [Pa] \n',...
        ref_sig(iMember), min_sig_all(iMember), max_sig_all(iMember));
end
fprintf(' ============================================= \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
% --->
figure;
plot(x_plot * (10^3), y_plot(2,:), 'r--', 'LineWidth', 1.0);
hold on;
plot(x_plot * (10^3), y_plot(3,:), 'r--', 'LineWidth', 1.0);
% plot(x_plot * (10^3), y_plot(1,:), 'r-', 'LineWidth', 1.0);
plot(min_eps * (10^3), min_sig, 'r^',...
    'LineWidth', 1.0, 'MarkerSize',12, 'MarkerFaceColor','w');
plot(max_eps * (10^3), max_sig, 'rv',...
    'LineWidth', 1.0, 'MarkerSize',12, 'MarkerFaceColor','w');
plot(ref_eps * (10^3), ref_sig, 'bs',...
    'LineWidth', 1.0, 'MarkerSize',9, 'MarkerFaceColor','w');
plot(data_x(1,:) * (10^3), data_x(2,:), 'ko',...
    'MarkerSize',3, 'MarkerFaceColor','w');
hold on;
set(gca,'FontName','Times');
set(gca,'FontSize',14);
set(gcf,'renderer','painters');
xlabel('Strain ($10^{-3}$ m/m)', 'Interpreter', 'latex');
ylabel('Stress ($10^{6}$ Pa)', 'Interpreter', 'latex');

figure;
hold on;
for iMember=1:nm
    plot(iMember, min_sig_all(iMember), 'r^',...
        'LineWidth', 1.0, 'MarkerSize',6, 'MarkerFaceColor','w');
    plot(iMember, max_sig_all(iMember), 'rv',...
        'LineWidth', 1.0, 'MarkerSize',6, 'MarkerFaceColor','w');
    plot(iMember, ref_sig(iMember), 'bs',...
        'LineWidth', 1.0, 'MarkerSize',6, 'MarkerFaceColor','w');
end
xlim( [0.5, iMember+0.5] );
set(gca,'FontName','Times');
set(gca,'FontSize',12);
set(gcf,'renderer','painters');
xlabel('Member index', 'Interpreter', 'latex');
ylabel('Stress ($10^{6}$ Pa)', 'Interpreter', 'latex');
if Flag.save == 1
    print('-depsc2', '-painters', '10bar_stress_bound');
end
% <---
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Flag.save == 1
    save( 'prog2_truss_stress_analysis.mat');
    diary off;
end

