%% Analysis of repeated experiments

% Let's look at the various statistical indicators
% The different experiments are used to define the error bounds

% See Quad. 2.3, 4/11/2021 for comprehensive plan


function [param_ok_new] = optimise(filename, combo, white_noise, mult_noise)

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
    elseif mult_noise == 3
      load(fullfile(parentdir,'data_both/',filename)); % up to c=1.68 (after the transition)
    end
end

val2 = [1.9:-0.002:1.68] - 1.78;

sol = sol(5000:end,1:end,1:end); % cut out transients


%% Analysis
% Estimate summary statistics

% Var
variance_matrix = reshape( var(sol,1), size(sol,2), size(sol,3) );  % over various experiments

% AC
AC_matrix = zeros(size(sol,2),size(sol,3));
for n=1:size(sol,3)
    for m=1:size(sol,2)
            coeff = corrcoef(sol(2:end,m,n) , sol(1:end-1,m,n));
            AC_matrix(m,n) = abs(coeff(1,2));
    end
end


% Coefficient Variation
mean_matrix = reshape( mean(sol,1), size(sol,2), size(sol,3) );
mean_mean = mean(mean_matrix')';
std_mean = std(mean_matrix');

% Shannon entropy
entropy_matrix1 = zeros(size(sol,2),size(sol,3)); % classical Shannon entropy
for n=1:size(sol,2)
    for m=1:size(sol,3)
        [h1,edges] = histcounts(sol(:,n,m),'Normalization', 'Probability'); % for Shannon entropy
        h1 = h1((h1~=0));
        entropy_matrix1(n,m) = -sum(h1.*log2(h1))*(edges(2)-edges(1));
    end
end

%%
param_ok_new = zeros(1,size(combo,2));
for index=1:size(combo,2)
    new_ind = combo(1,index) .*  variance_matrix + combo(2,index) .* AC_matrix + combo(3,index) * entropy_matrix1;
    param_ok_new(index) = testsignificance(new_ind,val2);
end

%%
% Let's have a look inside

% index = [70,110];
% new_ind = combo(1,index(1)) .*  variance_matrix + combo(2,index(1)) .* AC_matrix + combo(3,index(1)) * entropy_matrix1;
% mean_new_ind = mean(new_ind,2);
% std_mean = std(new_ind');
% new_ind2 = combo(1,index(2)) .*  variance_matrix + combo(2,index(2)) .* AC_matrix +  combo(3,index(2)) * entropy_matrix1;
% mean_new_ind2 = mean(new_ind2,2);
% std_mean2 = std(new_ind2');
% 
% figure()
% hold on
% plot(val2,mean_new_ind)
% plot(val2,mean_new_ind2)

end