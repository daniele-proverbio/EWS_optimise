%% Loop the analysis for all files, corresponding to different noise levels


clear all; close all; clc;

% Load files to analyse
parentdir =  fileparts(pwd);

white_noise = false;  % Do we want to analyse simulations with  white or multiplicative noise?
mult_noise = 2;       % 0 = sigma x; 1 = sigma x^2; 2 = sigma x^2/(1+x^2)

if white_noise
    files = dir(fullfile(parentdir,'data/*.mat'));   %files to analyse    
else
    if mult_noise == 0
        files = dir(fullfile(parentdir,'data_mn/*.mat'));   %files to analyse
    elseif mult_noise == 1
        files = dir(fullfile(parentdir,'data_mn1/*.mat'));   %files to analyse
    elseif mult_noise == 2
        files = dir(fullfile(parentdir,'data_mn2/*.mat'));   %files to analyse
    end
end

%%
params_ok = zeros(size(files,1),7);  % parameter values for which the increase is significant, for 7 statistical indicators [var,ac1,skewness,kurtosis,CV,index dispersion, shannon entropy]
counter = zeros(size(files,1),60);   % how many systems tipped at certain parameter values, for all considered parameter values from baseline to bifurcation point
val2 = [1.9:-0.002:1.68] - 1.78;     % parameter values considered (before and after transition)
checkplot = 0;                       % shall I plot an example figure?
sigmas = [0.01,0.012,0.014,0.016,0.018,0.02,0.025,0.03,0.035,0.04,0.045,0.05]; % noise level tested

% Call analysis function, looping on all files

for n=1:size(files,1)
    if n==2         % plot example figure only for sigma = 0.02
        checkplot = 1;
    else
        checkplot = 0;
    end
    [params_ok(n,:), counter(n,:)] = analysis(files(n).name,checkplot, white_noise, mult_noise);  % file name; shall I make the plots [1=yes, 0=no]?
end


%%  Plot params_ok
% For each statistical indicator

figure()
plot(sigmas,params_ok,'-o',linewidth=1.5)
legend({'Var','AC1','Skew','Kurt','CV','ID','H_S'},FontSize=20,Location='southeast',NumColumns=3)
ax = gca;
ax.FontSize = 20; 
ylabel('$c_{sig}$',fontsize=36,Interpreter='latex')
xlabel('$\sigma$',fontsize=36,Interpreter='latex')

%% Plot counters

figure()
imagesc(val2(1:60),sigmas,counter)
set(gca,'YDir','normal')
ax = gca;
ax.FontSize = 20; 
xlabel('$c-c_0$',fontsize=36,Interpreter='latex')
ylabel('$\sigma$',fontsize=36,Interpreter='latex')
title('Counts [%]')
colorbar


%% Figure example

analysis(files(2).name,1,white_noise,mult_noise);

