close all;
clear;
%
param_scale = 100;
domain_x   = 0.5*10^(-3) * param_scale;
num.sample = 100;
const.Young   = 20.0;
const.Poisson = 0.3;
const.std_Young   = 0.01 * const.Young;
const.std_Poisson = 0.01 * const.Poisson;


rng(33,'twister');
w = domain_x;
list_eps = -w + ( 2 * w * rand(3, num.sample) );
list_eps(3,:) = 2 * list_eps(3,:);
[~,idx_cur] = sort(list_eps(1,:));
list_eps = list_eps(:,idx_cur);
clear w idx_cur

rng(44,'twister');
list_Young = const.Young + (const.std_Young .* randn(num.sample, 1));
rng(55,'twister');
list_Poisson = const.Poisson + (const.std_Poisson .* randn(num.sample, 1));

list_sigma = zeros(3, num.sample);
for i=1:num.sample
    eps_cur = list_eps(:,i);
    dee = (  list_Young(i) / (1 - (list_Poisson(i)^2))  ) *...
        [1, list_Poisson(i), 0;...
        list_Poisson(i), 1, 0;...
        0, 0, (1/2) * (1 - list_Poisson(i))];
    list_sigma(:,i) = dee * eps_cur;
end

mean_Young   = mean(list_Young);
mean_Poisson = mean(list_Poisson);


fprintf(' Young modulus \n' );
fprintf('   mean = %2.5e [GPa] ; std_variation = %2.5e [GPa] \n',...
    mean(list_Young) * 10^(-1), std(list_Young) * 10^(-1));
fprintf(' Poisson ratio \n' );
fprintf('   mean = %2.5e ; std_variation = %2.5e \n',...
    mean(list_Poisson), std(list_Poisson));



list_norm_d_eps = zeros(1, num.sample);
list_norm_d_sig = zeros(1, num.sample);
for i=1:num.sample
    list_norm_d_eps(i) = mk_norm_dev_strain(list_eps(:,i), list_Poisson(i));
    list_norm_d_sig(i) = mk_norm_dev_sigma(list_sigma(:,i));
end





plot3(list_eps(1,:) / param_scale,...
    list_eps(2,:) / param_scale,...
    (1/2)*list_eps(3,:) / param_scale,...
    'ko', 'MarkerFaceColor','w', 'MarkerSize',3);
hold on;
grid on;
axis equal;
xlabel('$\varepsilon_{11}$ (m/m)', 'Interpreter', 'latex');
ylabel('$\varepsilon_{22}$ (m/m)', 'Interpreter', 'latex');
zlabel('$\varepsilon_{12}$ (m/m)', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',14);
print('-depsc2', '-painters', 'lin_isotropic_data_set_strain');
% savefig('lin_isotropic_data_set_strain.fig');


figure;
plot3(list_sigma(1,:), list_sigma(2,:), list_sigma(3,:), 'ko',...
    'MarkerFaceColor','w', 'MarkerSize',3);
hold on;
grid on;
axis equal;
xlabel('$\sigma_{11}$ (MPa)', 'Interpreter', 'latex');
ylabel('$\sigma_{22}$ (MPa)', 'Interpreter', 'latex');
zlabel('$\sigma_{12}$ (MPa)', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',14);
print('-depsc2', '-painters', 'lin_isotropic_data_set_stress');
% savefig('lin_isotropic_data_set_stress.fig');



figure;
plot( (1:num.sample), list_Young * param_scale, 'b.');
hold on;
xlabel('Sample ID', 'Interpreter', 'latex');
ylabel('Young''s modulus (MPa)', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
print('-depsc2', '-painters', 'lin_isotropic_data_set_Young');


figure;
plot( (1:num.sample), list_Poisson, 'b.');
hold on;
xlabel('Sample ID', 'Interpreter', 'latex');
ylabel('Poisson''s ratio', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
print('-depsc2', '-painters', 'lin_isotropic_data_set_Poisson');



figure;
plot(list_norm_d_eps / param_scale, list_norm_d_sig, 'ko',...
    'MarkerFaceColor','w', 'MarkerSize',2);
hold on;
xlim([0, 1.1*max(list_norm_d_eps / param_scale)]);
ylim([0, 1.1*max(list_norm_d_sig)]);
xlabel('equivalent strain (m/m)', 'Interpreter', 'latex');
ylabel('equivalent stress (MPa)', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
print('-depsc2', '-painters', 'lin_isotropic_data_set_eq_stress_strain');


save('lin_isotropic_data_set.mat',...
    'list_eps', 'list_sigma', 'param_scale',...
    'mean_Young', 'mean_Poisson',...
    'list_norm_d_eps', 'list_norm_d_sig');

