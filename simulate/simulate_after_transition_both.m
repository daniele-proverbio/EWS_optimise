%% Simulate system up to after the transition
% Both types of noise together

%% Optimise EWS for detection of Fold CT
% Looking for te best combination of summary statistis to detect critical
% transitions driven by fold bifurcations, in a biological toy model
% (autoactivating feedback loop motif for gene regulation)

% See Quad 2.3, 4/8/2021 for explanations and reasoning
% See Quad. 2.3, 4/11/2021 for comprehensive plan

%% Author
% Daniele Proverbio, 04/08/2021
% daniele.proverbio@outlook.com / @uni.lu
% University of Luxembourg


%% Prepare env
clear; close all; clc;

%% Initialize
% Choose type of multiplicative noise
alpha = 0.9;  % Proportion of white noise vs multiplicative noise (of type sigma*x)

N_Exp = 200;          %repeated experiments

time_start = 1;     % Start time of simulations 
N= 10000;            % Time points
dt = 0.01;           % Time step 
tic = 0;

val1 = 0.1;               %0.1 normal % Basal expression (constant, within accepted range), normally = 0.1 ; otherwise = 0  
val2 = 1.9:-0.002:1.68;     % Max Production (control parameter)  (up to after the bifurcation)
noise = 0.05;             % Noise level (diffusion term) -> lower than basin height

dist = val2 - 1.8;

sol = zeros((N*dt+1)/dt,length(val2),N_Exp);

%% SDE simulator
% Euler Maruyama scheme
for experiment =1 : N_Exp                        %repeated experiments
    
    x_in = 1.3; %+ (2.25-1.3)*rand;             % Initial condition (on upper branch) .  (initial conditions: x_in = [2.25, 1.3])
    
    for m = 1 : length(val2)   
        tic = 0;        
        for p = time_start-1+dt : dt : N*dt+1
            tic = tic + 1;          
            if tic == 1
                sol(tic,m,experiment) = x_in;
            else
                f = protein_production(p-dt,sol(tic-1,m,experiment),val1,val2(m));
                sol(tic,m,experiment) = sol(tic-1,m,experiment) + f * dt + alpha*(noise*sqrt(dt)*randn) + (1-alpha)*(sol(tic-1,m,experiment)*noise*sqrt(dt)*randn + 0.5*(sol(tic-1,m,experiment)*noise^2)*((sqrt(dt)*randn)^2 - dt));  % White as well as multiplicative noise. Milstein method
            end
        end
    end
end


%save('exp_after_transition.mat','sol', '-v7.3')  this only had 1
%experiment, to start trying
parentdir =  fileparts(pwd);

save(fullfile(parentdir,'data_both/multiple_exps_after_ct_both_050_a09.mat'),'sol', '-v7.3')


%% Equation
% Simulate simple equation for autocatalytic protein production (deterministic part) (Sharma,
% 2015)

function dxdt = protein_production(t,x,val1,val2)

K = val1;    
c = val2;    

dxdt = K + (c*x*x)/(1+x*x) - x;

end
