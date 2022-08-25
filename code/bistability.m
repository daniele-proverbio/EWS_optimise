%% Hysteresis
% Plot evolution graph for forward and backward parameters to get a
% hysteresis plot.
% Data are saved and come from simulations of "../simulate_SDE/varying
% parameter.m"

clear all
close all
clc


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
                cmap = [cmap; [0,0,0]];             % If Stable point, color it in black
            else
                cmap = [cmap; [1,0,0]];             % If ustable, color it in red
            end
        end
    end
end


%% Figure bistability

figure;
set(gca,'FontSize',15)
hold on
scatter(c_vector,solutions,10,cmap,'filled');
xlabel("c",fontsize=22)
ylabel("x",fontsize=22)
xlim([1.2,3])

%% Figure I/O
x = 0:0.01:1.6;
y =  (x.^2)./(1+x.^2);
y1 =  1.77*(x.^2)./(1+x.^2);
y2 =  2*(x.^2)./(1+x.^2);
%y3 =  2.4*(x.^2)./(1+x.^2);

fig = figure;
left_color = [0, 0.4470, 0.7410];
right_color = [0.8500, 0.3250, 0.0980];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
set(gca,'FontSize',15)

yyaxis right
hold on 
plot(x,y,'-.',linewidth=1.5,color=[0.8500, 0.3250, 0.0980])
plot(x,y1,'--',linewidth=1.5,color=[0.8500, 0.3250, 0.0980])
plot(x,y2,'-',linewidth=1.5,color=[0.8500, 0.3250, 0.0980])
%plot(x,y3,'--',linewidth=1.2,color=[0.8500, 0.3250, 0.0980])
ylabel('Feedback function',fontsize=22)
xlabel('x',fontsize=22)
ylim([-0.1,1.6])

yyaxis left
plot(x, x - 0.1,'k',linewidth=1.5,color=[0, 0.4470, 0.7410])
ylim([-0.1,1.6])
ylabel('Feedforward function',fontsize=22)

%% Sigmoidal

x1 = 0:0.01:2;
y5 = x1.^5./(1+x1.^5);
plot(y5,'k',linewidth=4)



