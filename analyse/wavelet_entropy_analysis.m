%% Quick analysis on wavelet entropy 
% What it does right before and after the transition
% (to make sure it reaches a maximum or at least stabilises after the
% transition)

clc; clear all; close all;


%% Analysis

mult_noise = false;   % Do we want to analyse simulations with  white or multiplicative noise?

parentdir =  fileparts(pwd);

if mult_noise
files = dir(fullfile(parentdir,'data_mn/*.mat'));   %files to analyse
else
files = dir(fullfile(parentdir,'data/*.mat'));   %files to analyse 
end

params_ok = zeros(size(files,1),7);  % parameter values for which the increase is significant, for 7 statistical indicators [var,ac1,skewness,kurtosis,CV,index dispersion, shannon entropy]
sigmas = [0.01,0.012,0.014,0.016,0.018,0.02,0.025,0.03,0.035,0.04,0.045,0.05]; % noise level tested

%% Call analysis function, looping on all files
for n=1:size(files,1)
    [params_ok(n,:)] = wavelet_ent(files(n).name,mult_noise);  % file name; shall I make the plots [1=yes, 0=no]?
end


%%  Plot params_ok
% For detrended wavelet entropy

figure()
plot(sigmas,params_ok,'-o',linewidth=1.1)
legend({'E_w'},FontSize=14,Location='southeast',NumColumns=3)
ax = gca;
ax.FontSize = 16; 
ylabel('$c_{sig}$',fontsize=26,Interpreter='latex')
xlabel('$\sigma$',fontsize=26,Interpreter='latex')



%%



%function param_ok_wentr = wavelet_ent(filename, mult_noise)
%% Load data

parentdir =  fileparts(pwd);

if mult_noise
    load(fullfile(parentdir,'data_mn/',filename)); % up to c=1.68 (after the transition)
else
    load(fullfile(parentdir,'data/',files(1).name)); % up to c=1.68 (after the transition)
end

val2 = [1.9:-0.002:1.68] - 1.78;

%% Wavelet Entropy: estimate

entropy_matrix = zeros(size(sol,2),size(sol,3)); % entropy_subband from wavelet decomposition (wentropy MATLAB function, see readme for further details)
for n=1:size(sol,2)
    for m=1:size(sol,3)
        entropy_matrix(n,m) = wentropy(sol(:,n,m),'shannon'); % wavelet
    end
end
% wavelet
mean_entropy_w = mean(entropy_matrix')' ;
std_entropy_w = std(entropy_matrix');

% %% Plot sample
% 
% fig_cut = 1; % Just to pack the figure a bit
% fig_end = length(val2);
% sample = 41;  % One sample trajectory
% num_fig = 8; 
% 
% figure()
% high_i = mean_entropy_w+std_entropy_w';
% low_i = mean_entropy_w-std_entropy_w';
% hold on
% r = patch([val2(fig_cut:fig_end) fliplr(val2(fig_cut:fig_end))], [high_i(fig_cut:fig_end)' fliplr(low_i(fig_cut:fig_end)')], 'b','LineStyle','none');
% %plot(val2(fig_cut:fig_end),ID_matrix(fig_cut:fig_end, sample),'b')
% plot(val2(fig_cut:fig_end),mean_entropy_w(fig_cut:fig_end),'k','LineWidth',1)
% xline(0,'--','linewidth',1)
% alpha(r,.05)
% ax = gca;
% ax.FontSize = 14;
% xlabel("c - c_0",'fontsize',18)
% ylabel("E_w",'fontsize',18)

%% Detrend linear pattern

P = polyfit(val2(1:47),mean_entropy_w(1:47),1);
yfit = P(1)*val2(1:47)+P(2);
%figure()
%hold on;
%scatter(val2(1:47),mean_entropy_w(1:47));
%plot(val2(1:47),yfit)

yfit2 =  P(1)*val2(1:61)+P(2);
% figure()
% hold on;
% scatter(val2(1:61),mean_entropy_w(1:61));
% plot(val2(1:61),yfit2)

detrended_matrix = entropy_matrix(1:61,:)' -  yfit2;
mean_detrended_w = mean(detrended_matrix)' ;
std_detrended_w = std(detrended_matrix);

%% Sample detrended figure

figure()
high_i = mean_detrended_w+std_detrended_w';
low_i = mean_detrended_w-std_detrended_w';
hold on
r = patch([val2(1:61) fliplr(val2(1:61))], [high_i(1:61)' fliplr(low_i(1:61)')], 'b','LineStyle','none');
%plot(val2(fig_cut:fig_end),ID_matrix(fig_cut:fig_end, sample),'b')
plot(val2(1:61),mean_detrended_w,'k','LineWidth',1)
xline(0,'--','linewidth',1)
alpha(r,.05)
ax = gca;
ax.FontSize = 14;
xlabel("c - c_0",'fontsize',18)
ylabel("E_w, detrended",'fontsize',18)

%%
det_mat_T = detrended_matrix';

[~, pv] = ttest(det_mat_T(size(det_mat_T,1)-10,:) , det_mat_T(size(det_mat_T,1)-1,:));


%% Get significant parameter
param_ok_wentr = testsignificance(detrended_matrix',val2(1:61));

%end

