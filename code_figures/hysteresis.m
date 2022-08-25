%% Hysteresis
% Plot evolution graph for forward and backward parameters to get a
% hysteresis plot.
% Data are saved and come from simulations of "varying_parameter_NF.m.m"

clear all
close all
clc

%% Load calculated sims

param_f = load('parameters_forward.mat');
param_b = load('parameters_backward.mat');
sol_f = load('sol_forward.mat');
sol_b = load('sol_backward.mat');


%% Equilibria

syms x            % Working with symbolic manipulation (more precise)
K = 0.1;          % Basal expression (constant, within accepted range); normally, K = 0.1
c = 0:0.01:4.5;   % Max Production (control parameter)

solutions=[];   % Vector of equilibria solutions
c_vector=[];    % Vector of accepted control parameters
cmap = [];      % Vector of colors


for m = 1:length(c)     % Check all c
    
    f = K + (c(m)*(x)*(x))/(1+(x)*(x)) - (x);     % My equation (vector field)
    soly = vpasolve(f == 0, x);         % Look for roots
    f_prime = diff(f);                  % Estimate derivatives
    
    for n = 1:length(soly)
        if isreal(soly(n))              % Only real roots allowed, obviously
            solutions = [solutions,soly(n)];
            c_vector = [c_vector,c(m)];
            if(vpa(subs(f_prime,x,soly(n))) < 0 )   % Check for linear stability -> eigenvalue <> 0
                cmap = [cmap; [0,0,1]];             % If Stable point, color it in blue
            else
                cmap = [cmap; [1,0,0]];             % If ustable, color it in red
            end
        end
    end
end


%% Figure
figure;
hold on
scatter(c_vector,solutions,10,cmap,'filled');
plot(param_f.parameters_spanned,sol_f.sol, 'k', linewidth=2);
plot(param_b.parameters_spanned,sol_b.sol, linewidth=2);
xlabel("p",fontsize=22)
ylabel("x",fontsize=22)
xlim([1.2,3])
annotation('arrow',[0.3 0.4],[0.7   0.7]);
annotation('arrow',[0.4 0.3],[0.6   0.6],color=[0.9290, 0.6940, 0.1250] 	);
