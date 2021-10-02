close all;
clear;
%
param_scale = 100;
domain_x   = 1.0*10^(-3) * param_scale;
num.sample = 500;
const.Young   = 25.0;
const.Poisson = 0.3;
const.std_Young   = 0.005 * const.Young;
const.std_Poisson = 0.005 * const.Poisson;


rng(33,'twister');
w = domain_x;
list_eps = -w + ( 2 * w * rand(3, num.sample) );
list_eps(3,:) = 2 * list_eps(3,:);
[~,idx_cur] = sort(list_eps(1,:));
list_eps = list_eps(:,idx_cur);
clear w idx_cur

rng(321,'twister');
list_Young1 = const.Young + (const.std_Young .* randn(num.sample, 1));
mean_Young(1)   = mean(list_Young1);
rng(432,'twister');
list_Young2 = 0.20 *...
    (  const.Young + (const.std_Young .* randn(num.sample, 1))  );
mean_Young(2)   = mean(list_Young2);

rng(543,'twister');
list_Poisson1 = const.Poisson + (const.std_Poisson .* randn(num.sample, 1));
mean_Poisson(1) = mean(list_Poisson1);
rng(654,'twister');
list_Poisson2 = 0.8 *...
    const.Poisson + (const.std_Poisson .* randn(num.sample, 1));
mean_Poisson(2) = mean(list_Poisson2);



list_sigma = zeros(3, num.sample);
for i=1:num.sample
    eps_cur = list_eps(:,i);
    dee =...
        [list_Young1(i) / (1 - (list_Poisson1(i) * list_Poisson2(i))),...
        (list_Poisson1(i) * list_Young2(i)) / (1 - (list_Poisson1(i) * list_Poisson2(i))),...
        0;...
        (list_Poisson1(i) * list_Young2(i)) / (1 - (list_Poisson1(i) * list_Poisson2(i))),...
        list_Young2(i) / (1 - (list_Poisson1(i) * list_Poisson2(i))),...
        0;...
        0,...
        0,...
        (list_Young1(i) * list_Young2(i)) /...
            (  list_Young1(i) + list_Young2(i) + (2 * list_Poisson1(i) * list_Young2(i))  )];
    list_sigma(:,i) = dee * eps_cur;
end


list_norm_d_eps = zeros(1, num.sample);
list_norm_d_sig = zeros(1, num.sample);
for i=1:num.sample
    list_norm_d_eps(i) =...
        mk_norm_dev_strain(list_eps(:,i), list_Poisson1(i), list_Poisson2(i));
    list_norm_d_sig(i) = mk_norm_dev_sigma(list_sigma(:,i));
end


fprintf(' Young modulus \n' );
fprintf('   mean = %2.5e [GPa] ; std_variation = %2.5e [GPa] \n',...
    mean(list_Young1) * 10^(-1), std(list_Young1) * 10^(-1));
fprintf('   mean = %2.5e [GPa] ; std_variation = %2.5e [GPa] \n',...
    mean(list_Young2) * 10^(-1), std(list_Young2) * 10^(-1));
fprintf(' Poisson ratio \n' );
fprintf('   mean = %2.5e ; std_variation = %2.5e \n',...
    mean(list_Poisson1), std(list_Poisson1));
fprintf('   mean = %2.5e ; std_variation = %2.5e \n',...
    mean(list_Poisson2), std(list_Poisson2));



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
print('-depsc2', '-painters', 'prog6_lin_orthotropic_data_set_strain');


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
print('-depsc2', '-painters', 'prog6_lin_orthotropic_data_set_stress');


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
print('-depsc2', '-painters', 'prog6_lin_orthotropic_data_set_eq_stress_strain');


save('lin_orthotropic_data_set.mat',...
    'list_eps', 'list_sigma', 'param_scale',...
    'mean_Young', 'mean_Poisson',...
    'list_norm_d_eps', 'list_norm_d_sig');

