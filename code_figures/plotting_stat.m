%% Plotting 
% I only use variance as an example

clear 
close all
clc

%% Loading
% Load the statistics indicators calculated in analysis.m of the
% "simulations cooperativity" study

n2 = load('statistics2.mat');
n2_full = load('multiple_exps_ct_wn_n2');

%%
cc0 = flip(0:0.005:0.2);
stop=26;  % To avoid getting exactly into n-tipping only

%% Sigma (2 sigma, in this case)

figure()
hold on
plot(cc0(1:stop),n2.statistics.mean_var(1:stop),'k',linewidth=1.2)
yline(n2.statistics.mean_var(1) + 2 * n2.statistics.std_var(1),'--k',linewidth=1.2)
ax = gca;
ax.FontSize = 14;
ylabel('Var',fontsize=20,Interpreter='latex',FontWeight='bold')
xlabel('$c-c_0$',fontsize=20,Interpreter='latex',FontWeight='bold')

%% Kendall's tau

figure()
hold on
plot(cc0(1:stop),n2.statistics.mean_var(1:stop),'k',linewidth=1.2)
ax = gca;
ax.FontSize = 14;
ylabel('Var',fontsize=20,Interpreter='latex',FontWeight='bold')
xlabel('$c-c_0$',fontsize=20,Interpreter='latex',FontWeight='bold')
rectangle(LineStyle ='--',Position=[0.1 0.0012 0.04 0.0018])
rectangle(LineStyle ='-.',Position=[0.105 0.00103 0.04 0.0018])

kendall1 = corr(n2.statistics.mean_var(12:20),n2.statistics.mean_var(11:19),'type','Kendall');

txt = ['\tau = ' num2str(round(kendall1,2))];
text(0.15,0.0052,txt,FontSize=20)


%% p-value var

% To show the distributions around Variance average values, for two example
% c values
en.sol=n2_full.sol(5000:end,1:end,1:end);

% Var
variance_matrix = reshape( var(en.sol,1), size(en.sol,2), size(en.sol,3) );  % over various experiments
mean_variance = mean(variance_matrix,2);
std_variance = std(variance_matrix');

[h, pv] = ttest(variance_matrix(22,:) , variance_matrix(1,:));
p_value = abs(pv);

figure()
hold on
plot(cc0(1:stop),n2.statistics.mean_var(1:stop),'k',linewidth=1)
boxplot([variance_matrix(22,:)',variance_matrix(1,:)'], 'symbol', '',Positions=[cc0(22) cc0(1)])  % NB: outliers removed for visualization purposes
ylim([-0.2,5]*10^(-3))
ax = gca;
ax.FontSize = 14;
xticklabels({cc0(35),cc0(1)})
ylabel('Var',fontsize=20,Interpreter='latex',FontWeight='bold')
xlabel('$c-c_0$',fontsize=20,Interpreter='latex',FontWeight='bold')
set(findobj(gca,'type','line'),'linew',3)
txt1 = ['p-value = ' num2str(round(p_value,4))];
text(0.15,0.0042,txt1,FontSize=20)
