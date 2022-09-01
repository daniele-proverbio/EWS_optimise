%% Plot optimal weights
% Results are saved after each run of "full_optimise.m", for different
% alpha

alpha = 0:0.1:1;

w =  [0, 0, 1 ;
       0, 0, 1 ;
       0.1, 0.2, 0.7;
       0.4, 0.3, 0.3;
       0.1, 0.2, 0.6
       0.6, 0.4, 0;
       0.6, 0.3, 0.1;
       0.7, 0.2, 0.1;
       1,0,0
       0.7, 0.1, 0.2;
       0.9, 0.1, 0]';

w_bis = [0,0,1;
         0.1, 0, 0.9;
         0, 0.2, 0.8;
         0.2, 0.4, 0.4;
         0.2, 0.3, 0.5;
         0.6, 0.4, 0;
         0.4, 0.4, 0.2;
         0.7, 0.3, 0;
         1,0,0;
         0.7,0.1,0.2;
         0.9,0.1,0]';

figure()
hold on
plot(alpha,w_bis(1,:),'-o',linewidth=1.8,color=[0, 0.4470, 0.7410])
plot(alpha,w(1,:),'--o',linewidth=1.8,color=[0, 0.4470, 0.7410])
plot(alpha,w_bis(2,:),'-o',linewidth=1.8,color=[0.8500, 0.3250, 0.0980])
plot(alpha,w(2,:),'--o',linewidth=1.8,color=[0.8500, 0.3250, 0.0980])
plot(alpha,w_bis(3,:),'-o',linewidth=1.8,color=[0.9290, 0.6940, 0.1250])
plot(alpha,w(3,:),'--o',linewidth=1.8,color=[0.9290, 0.6940, 0.1250])
ax = gca;
ax.FontSize = 20; 
ylabel('$\mathbf{\hat{w}}$',fontsize=36,Interpreter='latex')
xlabel('$\alpha$',fontsize=36,Interpreter='latex')
legend('Var','','AC(1)','','H_S','Position',[0.754464285714286 0.419047619047619 0.136607142857143 0.151190476190476],fontsize=20) 
