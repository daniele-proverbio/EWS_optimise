%% Loop the analysis for all files, corresponding to different noise levels


clear all; close all; clc;

% Load files to analyse
parentdir =  fileparts(pwd);
files = dir(fullfile(parentdir,'data/*.mat'));   %files to analyse

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
    [params_ok(n,:), counter(n,:)] = analysis(files(n).name,checkplot);  % file name; shall I make the plots [1=yes, 0=no]?
end


%%  Plot params_ok
% For each statistical indicator

figure()
plot(sigmas,params_ok,'-o',linewidth=1.1)
legend({'Var','AC1','Skew','Kurt','CV','ID','H_S'},FontSize=14,Location='southeast',NumColumns=3)
ax = gca;
ax.FontSize = 16; 
ylabel('$c_{sig}$',fontsize=26,Interpreter='latex')
xlabel('$\sigma$',fontsize=26,Interpreter='latex')

%% Plot counters

figure()
imagesc(val2(1:60),sigmas,counter)
set(gca,'YDir','normal')
ax = gca;
ax.FontSize = 16; 
xlabel('$c-c_0$',fontsize=26,Interpreter='latex')
ylabel('$\sigma$',fontsize=26,Interpreter='latex')
title('Counts [%]')
colorbar


%% Figure example

analysis(files(n).name,checkplot);

