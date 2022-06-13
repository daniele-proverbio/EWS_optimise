%% Create plots
% From simulations obtained and saved from the other files


%% Prepare env
clear; close all; clc;


%% Plot bistable equilibria
% Create bifurcation diagrams for different values of the parameter k

load('matrix_equilibria_k05');
load('matrix_equilibria_k02');
load('matrix_equilibria_k01');

% Scatter plot
figure;
hold on
scatter(matrix_equilibria_K05(1,1:end),matrix_equilibria_K05(2,1:end),5,matrix_equilibria_K05(3:5,1:end)','filled')
scatter(matrix_equilibria_K02(1,1:end),matrix_equilibria_K02(2,1:end),5,matrix_equilibria_K02(3:5,1:end)','filled')
scatter(matrix_equilibria_K01(1,1:end),matrix_equilibria_K01(2,1:end),5,matrix_equilibria_K01(3:5,1:end)','filled')
title('Equilibria of the system','FontSize',15);
xlabel('Max production rate c','FontSize',14);
ylabel('TF-A Concentration at equilibrium','FontSize',14);


% Plot
% figure;
% hold on
% plot(matrix_equilibria_K05(1,1:end),matrix_equilibria_K05(2,1:end))
% plot(matrix_equilibria_K02(1,1:end),matrix_equilibria_K02(2,1:end))
% plot(matrix_equilibria_K01(1,1:end),matrix_equilibria_K01(2,1:end))
