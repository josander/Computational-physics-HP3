%% Homeproblem 3b
% Task 1
set(0, 'defaultTextInterpreter', 'latex'); 
clc
clear all

phi = dlmread('phi.data');

[x y] = meshgrid(0:1/(length(phi)-1):1,0:1/(length(phi)-1):1);
figure(1);
surf(x,y,phi','EdgeAlpha',0.02);
figure(2)
plot(0:1/(length(phi)-1):1,phi(:,(length(phi)-1)/2 +1))

%% Task 2
% Plot the obtained result and compare it to the result from E6
clc

e6 = dlmread('phi_E6.data');
xE6 = linspace(0, length(e6), length(e6));

figure(2);
clf
plot(xE6,e6);
xlim([0 length(e6)]);


%% Task 3
% Plot the exact solution and compare with the obtained result from the
% simulation

exact = dlmread('phi_exact_5000x5000.txt');
%phi = dlmread();
figure(3)
plot(linspace(0,1,length(exact)),exact, linspace(0,1,length(phi)),phi(:,(length(phi)-1)/2 +1));



%% Plots task 1
set(0, 'defaultTextInterpreter', 'latex'); 
figure(1)
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
phi = dlmread('phi.data');

plot(0:1/(length(phi)-1):1,phi(:,(length(phi)-1)/2 +1))
y = ylabel('y [V]', 'fontsize', 12);   
x = xlabel('x [m]', 'fontsize', 12);   

plotTickLatex2D
print(gcf,'-depsc2','twoGrid.eps')

