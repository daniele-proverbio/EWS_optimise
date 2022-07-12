
%% Analyse data from Dai et al., 2012
% Reproduce, re-analyse and reinterpret them in light of my analysis (Ch 5 of thesis)
% Details in Quad 2.4, 8/7/2022

clc; clear all; close all;

%%
path = "/Users/daniele.proverbio/Documents/PHD_Lussemburgo/PHD/DRC/Reports/manuscript1_dynamical_context/doi_10.5061_dryad.p2481134__v1/data_indicator/";
file_names = ["subset_dilution250.txt", "subset_dilution500.txt", "subset_dilution750.txt", "subset_dilution1000.txt", "subset_dilution1133.txt", "subset_dilution1266.txt", "subset_dilution1400.txt", "subset_dilution1600.txt"];

dilution = [250,500,750,1000,1133,1266,1400,1600];

how_many = 50; % boostrapping repetitions
mean_mat = zeros(how_many,8); std_mat = zeros(how_many,8); var_mat = zeros(how_many,8); skew_mat = zeros(how_many,8); kurt_mat = zeros(how_many,8); cv_mat = zeros(how_many,8); id_mat = zeros(how_many,8); au_mat = zeros(how_many,8); corr_time_mat = zeros(how_many,8); entropy_mat = zeros(how_many,8);

% Bootstrap them all!
for k = 1:8
    [mean_mat(:,k), std_mat(:,k), var_mat(:,k), skew_mat(:,k), kurt_mat(:,k), cv_mat(:,k), id_mat(:,k), au_mat(:,k), corr_time_mat(:,k), entropy_mat(:,k)] = bootstrap(path, file_names(k),how_many);
end

%% Reproduce statistical analysis
% I plot those considered in the original paper, to check for consistency

% NB: they plot the standard error (std/sqrt(n)) before grouping the days;
% I plot the standard deviation after grouping. The reason is, I want to
% have a better feeling of the underlying distribution. The results are
% nonetheless consistent

% The mean is useful to assess the equilibrium values 
mean_mean = mean(mean_mat);
mean_std = std(mean_mat);

% Calculate indicators

std_mean =  mean(std_mat);
std_std = std(std_mat);

var_mean = mean(var_mat);
var_std = std(var_mat);

skew_mean = mean(skew_mat);
skew_std = std(skew_mat);

kurt_mean = mean(kurt_mat);
kurt_std = std(skew_mat);

cv_mean = mean(cv_mat);
cv_std = std(cv_mat);

id_mean = mean(id_mat);
id_std = std(id_mat);

au_mean = mean(au_mat);
au_std = std(au_mat);

ct_mean = mean(corr_time_mat);
ct_std = std(corr_time_mat);

entropy_mean = mean(entropy_mat);
entropy_std = std(entropy_mat);

%% Test significance

pv_std = testsignificance(std_mat);
pv_var = testsignificance(var_mat);
pv_cv = testsignificance(cv_mat);
pv_id = testsignificance(id_mat);
pv_entropy = testsignificance(entropy_mat);


%% Plotting

subplot(3,3,1)
errorbar(dilution, std_mean, std_std,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('Standard deviation (cells/{\mu}l)', fontsize = 12)
ylim([0,0.2*10^5])
xlim([0,1800])
yyaxis right
semilogy(dilution(2:end),pv_std,'-s');
yline(0.01,'-')
ylim([10^(-40),100])
ylabel('p-value')

subplot(3,3,2)
errorbar(dilution, var_mean, var_std,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('Variance (cells/{\mu}l)^2', fontsize = 12)
%ylim([0,0.2*10^5])
xlim([0,1800])
yyaxis right
semilogy(dilution(2:end),pv_var,'-s');
yline(0.01,'-')
ylim([10^(-40),100])
ylabel('p-value')

subplot(3,3,7)
errorbar(dilution, skew_mean, skew_std,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('Skewness', fontsize = 12)
ylim([-1,1])
xlim([0,1800])

subplot(3,3,8)
errorbar(dilution, kurt_mean, kurt_std,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('Kurtosis', fontsize = 12)
%ylim([-1,1])
xlim([0,1800])

subplot(3,3,4)
errorbar(dilution, cv_mean, cv_std,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('Coefficient Variation', fontsize = 12)
ylim([0,0.25])
xlim([0,1800])
yyaxis right
semilogy(dilution(2:end),pv_cv,'-s');
yline(0.01,'-')
ylim([10^(-40),100])
ylabel('p-value')

subplot(3,3,5)
errorbar(dilution, id_mean, id_std,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('Index dispersion', fontsize = 12)
%ylim([-1,1])
xlim([0,1800])
yyaxis right
semilogy(dilution(2:end),pv_id,'-s');
yline(0.01,'-')
ylim([10^(-40),100])
ylabel('p-value')

subplot(3,3,9)
errorbar(dilution, au_mean, au_std,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('AC(1)', fontsize = 12)
ylim([0,1])
xlim([0,1800])

subplot(3,3,6)
errorbar(dilution, ct_mean, ct_std,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('Autocorrelation time (days)', fontsize = 12)
ylim([0,15])
xlim([0,1800])

subplot(3,3,3)
errorbar(dilution, entropy_mean, entropy_std,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('Entropy', fontsize = 12)
%ylim([-1,1])
xlim([0,1800])
yyaxis right
semilogy(dilution(2:end),pv_entropy,'-s');
yline(0.01,'-')
ylim([10^(-40),100])
ylabel('p-value')



%% Multivariate indicators

% CV and entropy (a posteriori)
mult_mat1 = entropy_mat/max(entropy_mat(:,1)) + cv_mat/max(cv_mat(:,1));
mult_mean1 = mean(mult_mat1);
mult_std1 = std(mult_mat1);

% Var and entropy (same weights)
mult_mat2 = entropy_mat/max(entropy_mat(:,1)) + var_mat/max(var_mat(:,1));
mult_mean2 = mean(mult_mat2);
mult_std2 = std(mult_mat2);

% Var and entropy (more weights on entropy)
mult_mat3 = 0.9*entropy_mat/max(entropy_mat(:,1)) + 0.1*var_mat/max(var_mat(:,1));
mult_mean3 = mean(mult_mat3);
mult_std3 = std(mult_mat3);

% p-values
pv_mult1 = testsignificance(mult_mat1);
pv_mult2 = testsignificance(mult_mat2);
pv_mult3 = testsignificance(mult_mat3);

subplot(1,3,1)
errorbar(dilution, mult_mean1, mult_std1,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('CV + Entropy', fontsize = 12)
%ylim([-1,1])
xlim([0,1800])
yyaxis right
semilogy(dilution(2:end),pv_mult1,'-s');
yline(0.01,'-')
ylim([10^(-40),100])
ylabel('p-value')

subplot(1,3,2)
errorbar(dilution, mult_mean2, mult_std2,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('Variance + Entropy', fontsize = 12)
%ylim([-1,1])
xlim([0,1800])
yyaxis right
semilogy(dilution(2:end),pv_mult2,'-s');
yline(0.01,'-')
ylim([10^(-40),100])
ylabel('p-value')

subplot(1,3,3)
errorbar(dilution, mult_mean3, mult_std3,'-o')
xlabel('Dilution factor', fontsize = 12)
ylabel('0.1 Variance + 0.9 Entropy', fontsize = 12)
%ylim([-1,1])
xlim([0,1800])
yyaxis right
semilogy(dilution(2:end),pv_mult3,'-s');
yline(0.01,'-')
ylim([10^(-40),100])
ylabel('p-value')


%% Calculate p-values
% reference: dil = 250 (fartest from tipping point)

function p_value = testsignificance(A)
    %[~, pv] = ttest2(A(:,1) , A(:,3),'Vartype','unequal')

    p_value = zeros(7,1);  % comparison with no change at all
    for n=1:7
                [~, pv] = ttest2(A(:,1) , A(:,n+1),'Vartype','unequal');
                p_value(n) = abs(pv);
    end
end
