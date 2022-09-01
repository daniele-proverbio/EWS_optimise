%% Full optimise
% Loop for evey sigma

clear all; close all; clc;

% Load files to analyse
parentdir =  fileparts(pwd);

white_noise = true;  % Do we want to analyse simulations with white or multiplicative noise?
mult_noise = 3;       % 0 = sigma x; 1 = sigma x^2; 2 = sigma x^2/(1+x^2); 3 = both mult and white

if white_noise
    files = dir(fullfile(parentdir,'data/*.mat'));   %files to analyse    
else
    if mult_noise == 0
        files = dir(fullfile(parentdir,'data_mn/*.mat'));   %files to analyse
    elseif mult_noise == 1
        files = dir(fullfile(parentdir,'data_mn1/*.mat'));   %files to analyse
    elseif mult_noise == 2
        files = dir(fullfile(parentdir,'data_mn2/*.mat'));   %files to analyse
    elseif mult_noise == 3
        files = dir(fullfile(parentdir,'data_both/*_a02.mat'));   %files to analyse
    end
end


 %% Create combo of weights and run analysis

a = 0:0.1:1;
b = combvec(a,a,a);
c = b(1,:) + b(2,:) + b(3,:);
b = b(:,(c==1));

params_optimise = [];

val2 = [1.9:-0.002:1.68] - 1.78;     % parameter values considered (before and after transition)
sigmas = [0.01,0.012,0.014,0.016,0.018,0.02,0.025,0.03,0.035,0.04,0.045,0.05]; % noise level tested

for n=1:size(files,1)
    params_optimise = [params_optimise ; optimise(files(n).name,b,white_noise,mult_noise)];
end 

%% Scores

scores = sum(params_optimise,1);

figure('Position', [10 10 800 300])
plot(scores,linewidth=1.5)
ax = gca;
ax.FontSize = 20; 
ylabel('scores $\mathrm{S}$',fontsize=36,Interpreter='latex')
xlabel('Combination $\#$',fontsize=36,Interpreter='latex')
% if white_noise == true
%     title("WN")
% else
%     title("MN")
% end

%% Combinations
x = 1:64;
y = 1:3;

figure('Position', [10 10 800 300])
imagesc(x,y,b)
%set(gca,'YDir','normal')
ax = gca;
ax.FontSize = 23; 
xlabel('Combination $\#$',fontsize=36,Interpreter='latex')
%ylabel('$\sigma$',fontsize=36,Interpreter='latex')
yticks([1 2 3])
yticklabels({'AC(1)','Var','H_S'})
ytickangle(90)
colorbar


%% Uncertainty for sub-optimal combinations
% Vary weights up to 20%; check corresponding S on graph; relative difference from maximum
% Repeat the same for non-optimal weights (in this case, look for minimum
% Do it for all
% Results on Quad 2.4, 4/7/2022
