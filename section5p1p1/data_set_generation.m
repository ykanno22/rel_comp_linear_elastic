close all
clear
%
domain_x = 6*10^(-3);
num.division = 59;

const.Young = 1.0 * (10^3);
const.std = 1.0*10^(-1);

list_of_x = (-domain_x:((2*domain_x)/num.division):domain_x)';
num.sample    = length(list_of_x);

rng(22,'twister');
w = domain_x;
list_of_noisy_x = -w + (2.0 * w * rand(num.division+1 ,1));
list_of_noisy_x = sort(list_of_noisy_x);
clear w
rng(55,'twister');
% list_Young = const.Young + (const.std * randn(num.sample,1));
r = const.std * randn(num.sample,1);

list_of_noisy_f = zeros(num.sample,1);
for i=1:num.sample
    list_of_noisy_f(i) = ( const.Young * list_of_noisy_x(i) ) + r(i);
end

samples_eps = list_of_noisy_x;
samples_sig = list_of_noisy_f;
list_Young  = samples_sig ./ samples_eps;
mean_Young  = mean(list_Young);

fprintf(' Young modulus \n' );
fprintf('   mean = %2.4e [GPa] ; std_variation = %2.4e [GPa] \n',...
    mean(list_Young) * 10^(-3), std(list_Young) * 10^(-3));


figure;
plot(samples_eps * (10^3), samples_sig, 'ko',...
    'MarkerFaceColor','w', 'MarkerSize',3)
hold on;
set(gca,'FontName','Times');
set(gca,'FontSize',14);
set(gcf,'renderer','painters');
ylim( [-7,7] );
xlabel('Strain ($10^{-3}$ m/m)', 'Interpreter', 'latex');
ylabel('Stress ($10^{6}$ Pa)', 'Interpreter', 'latex');
print('-depsc2', '-painters', 'prog2_truss_data_set');


save('data_set.mat', 'samples_eps', 'samples_sig', 'mean_Young');


