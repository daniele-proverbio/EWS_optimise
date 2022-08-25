%% Simulate SDE with a slowly varying parameter
% Plot how retrospective analysis might fail in distinguishing abrupt vs
% continuous transition

%% Author
% Daniele Proverbio, 14/01/2022
% daniele.proverbio@outlook.com
% University of Luxembourg

%% Prepare env
clear;  clc;
close all;

%% Initialize

time_start = 1;       % Start time of simulations 
time_stop = 1200;     % End recording simulations 
dt = 0.5;             % Time step
sol = zeros((time_stop-time_start),1);
parameters_spanned = zeros(time_stop-time_start,1);
tic = 0;

x_in1 = 2.16;            % Initial condition (on upper branch) .  
x_in2 = 1.1;            % Initial condition (on lower branch) . 
vals1 = [0.01, 0.5];     % vals1(1) = Basal expression (constant, within accepted range)    % vals1(2) = Max Production (control parameter)  % Moving from right to left 
vals2 = [0.01, -0.5];     % vals2(1) = Basal expression (constant, within accepted range)    % vals2(2) = Max Production (control parameter)  % Moving from left to right


noise = 0.02;           % Noise level (diffusion term) -> pretty low, at least lower than basin height

% SDE simulator
% Euler Maruyama scheme
for p = time_start:dt:time_stop
    tic = tic + 1;
    
    if tic == 1
        sol(tic) = x_in2;
        parameters_spanned(1) = vals2(2);
    else
        vals2(2) = vals2(2) + 0.001;     % Increase/descrease (+/-) the value of the parameter (depending on forward or backward trend)
        parameters_spanned(tic) = vals2(2);
        f = normal_form(p-dt,sol(tic-1),vals2);
        sol(tic) = sol(tic-1) + f * dt + noise*sqrt(dt)*randn;
    end
    
end


%% Bif diagram

% Perform calculation

syms x            % Working with symbolic manipulation (more precise)

c = -0.5:0.001:2;   % control parameter (corresponds to beta

solutions=[];   % Vector of equilibria solutions
c_vector=[];    % Vector of accepted control parameters
cmap = [];      % Vector of colors


for m = 1:length(c)     % Check all c
    
    f = vals2(1) + c(m)*(x-1) - (x-1)^3  ;      % My equation (vector field): substitute with desired Normal Form
    
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




%% Original figure
figure; 
hold on
%scatter(c_vector,solutions,3,cmap,'filled');
plot(parameters_spanned(:),sol,Color='k',linewidth=2);
xlabel("p",fontsize=22)
ylabel("x",fontsize=22)
hold off

%% Equation
% Supercritical pitchfork

function dxdt = normal_form(t,x,vals)   

K = vals(1);    
c = vals(2);    

dxdt = K + c*(x-1) - (x-1)^3;  % Slight translation, otherwise x=0 remains unchanged

end
