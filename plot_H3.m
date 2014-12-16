%% Homeproblem 3b

clc
clear all

phi = dlmread('phi.data');

x = linspace(0,1,length(phi));
y = linspace(0,1,length(phi));

figure(1);
plot3(x,y,phi);


%% Task 3
% Plot the exact solution and compare with the obtained result

exact = dlmread('exact.data');
phi = dlmread();

plot(exact);

