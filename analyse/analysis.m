%% Analysis of repeated experiments

% Let's look at the various statistical indicators
% The different experiments are used to define the error bounds

% See Quad. 2.3, 4/11/2021 for comprehensive plan


function [params_ok, counter] = analysis(filename, plotcheck, white_noise, mult_noise)

%% Prepare

parentdir =  fileparts(pwd);

if white_noise
      load(fullfile(parentdir,'data/',filename)); % up to c=1.68 (after the transition)
else
    if mult_noise == 0
      load(fullfile(parentdir,'data_mn/',filename)); % up to c=1.68 (after the transition)
    elseif mult_noise == 1
      load(fullfile(parentdir,'data_mn1/',filename)); % up to c=1.68 (after the transition)
    elseif mult_noise == 2
      load(fullfile(parentdir,'data_mn2/',filename)); % up to c=1.68 (after the transition)
    end
end

val2 = [1.9:-0.002:1.68] - 1.78;

sol = sol(5000:end,1:end,1:end); % cut out transients

%% Analysis 
% Estimate summary statistics

% Var
variance_matrix = reshape( var(sol,1), size(sol,2), size(sol,3) );  % over various experiments
mean_variance = mean(variance_matrix,2);
std_variance = std(variance_matrix');

% AC
AC_matrix = zeros(size(sol,2),size(sol,3));
for n=1:size(sol,3)
    for m=1:size(sol,2)
            coeff = corrcoef(sol(2:end,m,n) , sol(1:end-1,m,n));
            AC_matrix(m,n) = abs(coeff(1,2));
    end
end
mean_AC = mean(AC_matrix,2);
std_AC = std(AC_matrix');

% Skewness
skew_matrix = reshape( skewness(sol,1), size(sol,2), size(sol,3) );  % over various experiments
mean_skew = mean(skew_matrix,2);
std_skew = std(skew_matrix');

% Kurtosis
kurt_matrix = reshape( kurtosis(sol,1), size(sol,2), size(sol,3) );  % over various experiments
mean_kurt= mean(kurt_matrix,2);
std_kurt = std(kurt_matrix');

% Coefficient Variation
mean_matrix = reshape( mean(sol,1), size(sol,2), size(sol,3) );
mean_mean = mean(mean_matrix')';
std_mean = std(mean_matrix');
mean_CV = sqrt(mean_variance)./mean_mean;
std_CV = sqrt( (std_mean./(2*sqrt(mean_variance').*mean_mean' )).^2 + (sqrt(mean_variance').*std_mean./(mean_mean'.^2)).^2); % error propagation

% Index of dispersion
mean_ID = mean_variance./mean_mean;
std_ID = sqrt( (std_variance./mean_mean').^2 + (mean_variance'.*std_mean./(mean_mean'.^2)).^2); % error propagation

% Entropy
entropy_matrix = zeros(size(sol,2),size(sol,3)); % entropy_subband from wavelet decomposition (wentropy MATLAB function, see readme for further details)
entropy_matrix1 = zeros(size(sol,2),size(sol,3)); % classical Shannon entropy
for n=1:size(sol,2)
    for m=1:size(sol,3)
        entropy_matrix(n,m) = wentropy(sol(:,n,m),'shannon'); % wavelet

        [h1,edges] = histcounts(sol(:,n,m),'Normalization', 'Probability'); % for Shannon entropy
        h1 = h1((h1~=0));
        entropy_matrix1(n,m) = -sum(h1.*log2(h1))*(edges(2)-edges(1));
    end
end
% wavelet
mean_entropy_w = mean(entropy_matrix')' ;
std_entropy_w = std(entropy_matrix');
%shannon
mean_entropy_s = mean(entropy_matrix1')' ;
std_entropy_s = std(entropy_matrix1');


%% Plot benchmark
% Only for one example



if plotcheck == 1

    fig_cut = 1; % Just to pack the figure a bit
    fig_end = length(val2);
    sample = 41;  % One sample trajectory
    num_fig = 8; 
    ax_font = 18;
    label_font = 22;
    
    figure('position',[10 10 900 900])
    subplot(ceil(num_fig/2),2,1)
    hold on
    high = mean_variance+std_variance';
    low = scaleStd(mean_variance, std_variance, variance_matrix, fig_cut);
    p = patch([val2(fig_cut:fig_end) fliplr(val2(fig_cut:fig_end))], [high(fig_cut:fig_end)' fliplr(low(fig_cut:fig_end)')], 'b','LineStyle','none');
    %plot(val2(fig_cut:end),variance_matrix(fig_cut:end, sample),'b')
    plot(val2(fig_cut:fig_end),mean_variance(fig_cut:fig_end),'k','LineWidth',1)
    xline(0,'--','linewidth',1)
    alpha(p,.05)
    ax = gca;
    ax.FontSize = ax_font; 
    xlabel("c - c_0",'fontsize',label_font)
    ylabel("Var.",'fontsize',label_font)
    
    subplot(ceil(num_fig/2),2,2)
    hold on
    high_a = mean_AC+std_AC';
    low_a = mean_AC-std_AC';
    q = patch([val2(fig_cut:fig_end) fliplr(val2(fig_cut:fig_end))], [high_a(fig_cut:fig_end)' fliplr(low_a(fig_cut:fig_end)')], 'b','LineStyle','none');
    %plot(val2(fig_cut:fig_end),AC_matrix(fig_cut:fig_end, sample),'b')
    plot(val2(fig_cut:fig_end),mean_AC(fig_cut:fig_end),'k','LineWidth',1)
    xline(0,'--','linewidth',1)
    alpha(q,.05)
    ax = gca;
    ax.FontSize = ax_font; 
    xlabel("c - c_0",'fontsize',label_font)
    ylabel("AC(1)",'fontsize',label_font)
    
    subplot(ceil(num_fig/2),2,3)
    high_s = mean_skew+std_skew';
    low_s = mean_skew-std_skew';
    hold on
    r = patch([val2(fig_cut:fig_end) fliplr(val2(fig_cut:fig_end))], [high_s(fig_cut:fig_end)' fliplr(low_s(fig_cut:fig_end)')], 'b','LineStyle','none');
    %plot(val2(fig_cut:fig_end),skew_matrix(fig_cut:fig_end, sample),'b')
    plot(val2(fig_cut:fig_end),mean_skew(fig_cut:fig_end),'k','LineWidth',1)
    xline(0,'--','linewidth',1)
    alpha(r,.05)
    ax = gca;
    ax.FontSize = ax_font;
    xlabel("c - c_0",'fontsize',label_font)
    ylabel("Skew.",'fontsize',label_font)
    
    subplot(ceil(num_fig/2),2,4)
    high_s = mean_kurt+std_kurt';
    low_s = mean_kurt-std_kurt';
    hold on
    r = patch([val2(fig_cut:fig_end) fliplr(val2(fig_cut:fig_end))], [high_s(fig_cut:fig_end)' fliplr(low_s(fig_cut:fig_end)')], 'b','LineStyle','none');
    %plot(val2(fig_cut:fig_end),kurt_matrix(fig_cut:fig_end, sample),'b')
    plot(val2(fig_cut:fig_end),mean_kurt(fig_cut:fig_end),'k','LineWidth',1)
    xline(0,'--','linewidth',1)
    alpha(r,.05)
    ax = gca;
    ax.FontSize = ax_font;
    xlabel("c - c_0",'fontsize',label_font)
    ylabel("Kurt.",'fontsize',label_font)
    
    subplot(ceil(num_fig/2),2,5)
    high_c = mean_CV+std_CV';
    low_c = mean_CV-std_CV';
    hold on
    r = patch([val2(fig_cut:fig_end) fliplr(val2(fig_cut:fig_end))], [high_c(fig_cut:fig_end)' fliplr(low_c(fig_cut:fig_end)')], 'b','LineStyle','none');
    %plot(val2(fig_cut:fig_end),sqrt(variance_matrix(fig_cut:fig_end, sample))./mean_matrix(fig_cut:fig_end, sample),'b')
    plot(val2(fig_cut:fig_end),mean_CV(fig_cut:fig_end),'k','LineWidth',1)
    xline(0,'--','linewidth',1)
    alpha(r,.05)
    ax = gca;
    ax.FontSize = ax_font;
    xlabel("c - c_0",'fontsize',label_font)
    ylabel("CV",'fontsize',label_font)
    
    subplot(ceil(num_fig/2),2,6)
    high_i = mean_ID+std_ID';
    low_i = mean_ID-std_ID';
    hold on
    r = patch([val2(fig_cut:fig_end) fliplr(val2(fig_cut:fig_end))], [high_i(fig_cut:fig_end)' fliplr(low_i(fig_cut:fig_end)')], 'b','LineStyle','none');
    %plot(val2(fig_cut:fig_end),(variance_matrix(fig_cut:fig_end, sample))./mean_matrix(fig_cut:fig_end, sample),'b')
    plot(val2(fig_cut:fig_end),mean_ID(fig_cut:fig_end),'k','LineWidth',1)
    xline(0,'--','linewidth',1)
    alpha(r,.05)
    ax = gca;
    ax.FontSize = ax_font;
    xlabel("c - c_0",'fontsize',label_font)
    ylabel("ID",'fontsize',label_font)
    
%     subplot(ceil(num_fig/2),2,8)
%     high_i = mean_entropy_w+std_entropy_w';
%     low_i = mean_entropy_w-std_entropy_w';
%     hold on
%     r = patch([val2(fig_cut:fig_end) fliplr(val2(fig_cut:fig_end))], [high_i(fig_cut:fig_end)' fliplr(low_i(fig_cut:fig_end)')], 'b','LineStyle','none');
%     %plot(val2(fig_cut:fig_end),ID_matrix(fig_cut:fig_end, sample),'b')
%     plot(val2(fig_cut:fig_end),mean_entropy_w(fig_cut:fig_end),'k','LineWidth',1)
%     xline(0,'--','linewidth',1)
%     alpha(r,.05)
%     ax = gca;
%     ax.FontSize = ax_font;
%     xlabel("c - c_0",'fontsize',18)
%     ylabel("E_w",'fontsize',18)
    
    subplot(ceil(num_fig/2),2,7)
    high_i = mean_entropy_s+std_entropy_s';
    low_i = mean_entropy_s-std_entropy_s';
    hold on
    r = patch([val2(fig_cut:fig_end) fliplr(val2(fig_cut:fig_end))], [high_i(fig_cut:fig_end)' fliplr(low_i(fig_cut:fig_end)')], 'b','LineStyle','none');
    %plot(val2(fig_cut:fig_end),ID_matrix(fig_cut:fig_end, sample),'b')
    plot(val2(fig_cut:fig_end),mean_entropy_s(fig_cut:fig_end),'k','LineWidth',1)
    xline(0,'--','linewidth',1)
    alpha(r,.05)
    ax = gca;
    ax.FontSize = ax_font;
    xlabel("c - c_0",'fontsize',label_font)
    ylabel("H_s",'fontsize',label_font)

end

%% Are the changes significant?
% Analysis of p-values, and of when they pass a threshold p = 0.05.
% This estimates the lead "time"

param_ok_var = testsignificance(variance_matrix,val2);
param_ok_AC = testsignificance(AC_matrix,val2);
param_ok_skew = testsignificance(skew_matrix,val2);
param_ok_kurt = testsignificance(kurt_matrix,val2);
param_ok_cv = testsignificance(sqrt(variance_matrix)./mean_matrix,val2);
param_ok_id = testsignificance(variance_matrix./mean_matrix,val2);
param_ok_ent = testsignificance(entropy_matrix1,val2);

params_ok = [param_ok_var,param_ok_AC,param_ok_skew,param_ok_kurt,param_ok_cv,param_ok_id,param_ok_ent];

%% How many tippings occur at each param value
% Just up to bifurcation point
counter = zeros(60,1);

for n=1:60  % Just up to bif point (the value before; afetrwards, everything tips necessarily)
    for m=1:size(sol,3)
        if any(sol(:,n,m)<0.3)  % 0.3 is a good value to identify a strong transition between alternative states
            counter(n) = counter(n)+1;
        end
    end
end
counter = (counter/200)';  % normalise

end