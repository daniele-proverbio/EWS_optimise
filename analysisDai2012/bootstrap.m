%% Analyse data from Dai et al., 2012
% And reinterpret them in light of my analysis (Ch 5 of thesis)
% Details in Quad 2.4, 8/7/2022

function [mean_b, std_b, var_b, skew_b, kurt_b, cv_b, id_b, au_b, corr_time_b, entropy_b] = bootstrap(path,file,how_many)

A = readtable(path+file);
A = A{:,:};

cell_d = 1.1*10^5*A; % cell density

%% Bootstrapping procedure

rng('default')  % For reproducibility
mean_b = reshape(bootstrp(how_many/5,@mean,cell_d),[],1);       % Average (equilibrium)
std_b = reshape(bootstrp(how_many/5,@std,cell_d),[],1);         % Std
var_b = reshape(bootstrp(how_many/5,@var,cell_d),[],1);         % Var
skew_b = reshape(bootstrp(how_many/5,@skewness,cell_d),[],1);   % Skewness
kurt_b = reshape(bootstrp(how_many/5,@kurtosis,cell_d),[],1);   % Kurtosis
cv_b = std_b./mean_b;                                    % Coeff Var
id_b = var_b./mean_b;                                    % Index dispersion

% Custom replacement bootstrapping for custom functions
au_b = zeros(how_many,1);                         % AC1
corr_time_b = zeros(how_many,1);                  % Autocorrelation time
entropy_b = zeros(how_many,1);                    % Shannon entropy
for b=1:how_many
    k = ceil(rand(size(cell_d,1),1)*size(cell_d,1));
    mnew = cell_d(k,:);
    [au_b(b), corr_time_b(b)] = autocorr1(mnew);
    entropy_b(b) = entropy1(mnew);
end

end

%% Extra functions

function [autocorr1, corr_time] = autocorr1(A)
    autocorr1 = zeros(4,1);                      % AC1
    for a=1:4
        autocorr1(a) = corr(A(1:end,a), A(1:end,a+1));
        if autocorr1(a) < 0
            autocorr1(a) = 0;
        end
    end
    autocorr1 = mean(autocorr1);
    if autocorr1 > 0                             % correlation time (with cutoff, see SupMat formula 4)
        corr_time = -1/log(autocorr1);
    else 
        corr_time = 0;
    end
end

function entropy1 = entropy1(A)
    entropy1 = zeros(4,1);
    for a=1:4
        [h1,edges] = histcounts(A(1:end,a),'Normalization', 'Probability'); % for Shannon entropy
        h1 = h1((h1~=0));
        entropy1(a) = -sum(h1.*log2(h1))*(edges(2)-edges(1));
    end
    entropy1 = mean(entropy1);
end
