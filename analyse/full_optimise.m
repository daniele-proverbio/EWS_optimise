%% Full optimise
% Loop for evey sigma

clear all; close all; clc;

% Load files to analyse
parentdir =  fileparts(pwd);

mult_noise = true;   % Do we want to analyse simulations with  white or multiplicative noise?
 if mult_noise
    files = dir(fullfile(parentdir,'data_mn/*.mat'));   %files to analyse
 else
    files = dir(fullfile(parentdir,'data/*.mat'));   %files to analyse 
 end


 %% Create combo of weights and run analysis

a = 0:0.1:1;
b = combvec(a,a,a);
c = b(1,:) + b(2,:) + b(3,:) ;
b = b(:,(c==1));

params_optimise = [];

val2 = [1.9:-0.002:1.68] - 1.78;     % parameter values considered (before and after transition)
sigmas = [0.01,0.012,0.014,0.016,0.018,0.02,0.025,0.03,0.035,0.04,0.045,0.05]; % noise level tested

for n=1:size(files,1)
    params_optimise = [params_optimise ; optimise(files(n).name,b,mult_noise)];
end 

%%

scores = sum(params_optimise,1);

figure()
plot(scores,linewidth=1.2)
ax = gca;
ax.FontSize = 16; 
ylabel('scores $\mathrm{S}$',fontsize=26,Interpreter='latex')
xlabel('Combination $\#$',fontsize=26,Interpreter='latex')
if mult_noise == true
    title("MN")
else
    title("WN")
end
