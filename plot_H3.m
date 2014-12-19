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
xE6 = linspace(0, 1, length(e6));

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

%% Plot computational times
gridsizes = [21 41 81 161 321 641 1281];
compGS = [80142 857844 9941913 119528568];
compTG = [234123 6824452 215294115];
compMGV = [33219 113406 487866 2338847 9462117 38044587 152553457];
compMGW = [26404 113770 647508 2343230 9774510 39965070 161722190];
compFMG = [26520 91051 397394 1571256 6309938 25346240 101649762];
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);
loglog(gridsizes, compMGV,gridsizes, compMGW,gridsizes, compFMG, gridsizes(1:length(compGS)), compGS)
y = ylabel('Number of Gauss-Seidel iterations [ ]', 'fontsize', 12);   
x = xlabel('Gridsize J [ ]', 'fontsize', 12);   
plotTickLatex2D
l = legend('Multigrid V','Multigrid W','Full Multigrid','Gauss-Seidel');
set(l,'Location','East')
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);

print(gcf,'-depsc2','nComp.eps')


%% Plots task 1
set(0, 'defaultTextInterpreter', 'latex'); 
figure(1)
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);


plot(0:1/(length(phiTG21)-1):1,phi(:,(length(phiTG21)-1)/2 +1))
y = ylabel('$\Phi$(x,y) [V]', 'fontsize', 12);    
x = xlabel('x [m]', 'fontsize', 12);   

plotTickLatex2D

print(gcf,'-depsc2','twoGrid.eps')


%% Plots task 2
halfpoints = (gridsizes -1 )./2 +1;

%% Plot MGV
figure(2)
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

plot(0:1/(gridsizes(1)-1):1,phiMGV21(:,halfpoints(1)),0:1/(gridsizes(2)-1):1,phiMGV41(:,halfpoints(2)),0:1/(gridsizes(3)-1):1,phiMGV81(:,halfpoints(3)),0:1/(gridsizes(4)-1):1,phiMGV161(:,halfpoints(4)),0:1/(gridsizes(5)-1):1,phiMGV321(:,halfpoints(5)),0:1/(gridsizes(6)-1):1,phiMGV641(:,halfpoints(6)),0:1/(gridsizes(7)-1):1,phiMGV1281(:,halfpoints(7))) 
y = ylabel('$\Phi$(x,y) [V]', 'fontsize', 12);   
x = xlabel('x [m]', 'fontsize', 12);   
plotTickLatex2D
l = legend('J = 21','J = 41','J = 81','J = 161','J = 321','J = 641','J = 1281');
set(l,'Location','NorthWest')


set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
print(gcf,'-depsc2','MGV.eps')
clf
plot(0:1/(gridsizes(1)-1):1,phiMGV21(:,halfpoints(1)),0:1/(gridsizes(2)-1):1,phiMGV41(:,halfpoints(2)),0:1/(gridsizes(3)-1):1,phiMGV81(:,halfpoints(3)),0:1/(gridsizes(4)-1):1,phiMGV161(:,halfpoints(4)),0:1/(gridsizes(5)-1):1,phiMGV321(:,halfpoints(5)),0:1/(gridsizes(6)-1):1,phiMGV641(:,halfpoints(6)),0:1/(gridsizes(7)-1):1,phiMGV1281(:,halfpoints(7))) 
y = ylabel('$\Phi$(x,y) [V]', 'fontsize', 12);   
x = xlabel('x [m]', 'fontsize', 12);   
axis([0.5 0.7 0 1.25])
plotTickLatex2D

l = legend('J = 21','J = 41','J = 81','J = 161','J = 321','J = 641','J = 1281');


set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
print(gcf,'-depsc2','MGV2.eps')

%% Plot MGW
figure(2)
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

plot(0:1/(gridsizes(1)-1):1,phiMGW21(:,halfpoints(1)),0:1/(gridsizes(2)-1):1,phiMGW41(:,halfpoints(2)),0:1/(gridsizes(3)-1):1,phiMGW81(:,halfpoints(3)),0:1/(gridsizes(4)-1):1,phiMGW161(:,halfpoints(4)),0:1/(gridsizes(5)-1):1,phiMGW321(:,halfpoints(5)),0:1/(gridsizes(6)-1):1,phiMGW641(:,halfpoints(6)),0:1/(gridsizes(7)-1):1,phiMGW1281(:,halfpoints(7))) 
y = ylabel('$\Phi$(x,y) [V]', 'fontsize', 12);   
x = xlabel('x [m]', 'fontsize', 12);   
plotTickLatex2D
l = legend('J = 21','J = 41','J = 81','J = 161','J = 321','J = 641','J = 1281');
set(l,'Location','NorthWest')


set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
print(gcf,'-depsc2','MGW.eps')
clf
plot(0:1/(gridsizes(1)-1):1,phiMGW21(:,halfpoints(1)),0:1/(gridsizes(2)-1):1,phiMGW41(:,halfpoints(2)),0:1/(gridsizes(3)-1):1,phiMGW81(:,halfpoints(3)),0:1/(gridsizes(4)-1):1,phiMGW161(:,halfpoints(4)),0:1/(gridsizes(5)-1):1,phiMGW321(:,halfpoints(5)),0:1/(gridsizes(6)-1):1,phiMGW641(:,halfpoints(6)),0:1/(gridsizes(7)-1):1,phiMGW1281(:,halfpoints(7))) 
y = ylabel('$\Phi$(x,y) [V]', 'fontsize', 12);   
x = xlabel('x [m]', 'fontsize', 12);   
axis([0.5 0.7 0 1.25])
plotTickLatex2D

l = legend('J = 21','J = 41','J = 81','J = 161','J = 321','J = 641','J = 1281');


set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
print(gcf,'-depsc2','MGW2.eps')

%% Plot FMG
figure(2)
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

plot(0:1/(gridsizes(1)-1):1,phiFMG21(:,halfpoints(1)),0:1/(gridsizes(2)-1):1,phiFMG41(:,halfpoints(2)),0:1/(gridsizes(3)-1):1,phiFMG81(:,halfpoints(3)),0:1/(gridsizes(4)-1):1,phiFMG161(:,halfpoints(4)),0:1/(gridsizes(5)-1):1,phiFMG321(:,halfpoints(5)),0:1/(gridsizes(6)-1):1,phiFMG641(:,halfpoints(6)),0:1/(gridsizes(7)-1):1,phiFMG1281(:,halfpoints(7))) 
y = ylabel('$\Phi$(x,y) [V]', 'fontsize', 12);   
x = xlabel('x [m]', 'fontsize', 12);   
plotTickLatex2D
l = legend('J = 21','J = 41','J = 81','J = 161','J = 321','J = 641','J = 1281');
set(l,'Location','NorthWest')


set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
print(gcf,'-depsc2','FMG.eps')
clf
plot(0:1/(gridsizes(1)-1):1,phiFMG21(:,halfpoints(1)),0:1/(gridsizes(2)-1):1,phiFMG41(:,halfpoints(2)),0:1/(gridsizes(3)-1):1,phiFMG81(:,halfpoints(3)),0:1/(gridsizes(4)-1):1,phiFMG161(:,halfpoints(4)),0:1/(gridsizes(5)-1):1,phiFMG321(:,halfpoints(5)),0:1/(gridsizes(6)-1):1,phiFMG641(:,halfpoints(6)),0:1/(gridsizes(7)-1):1,phiFMG1281(:,halfpoints(7))) 
y = ylabel('$\Phi$(x,y) [V]', 'fontsize', 12);   
x = xlabel('x [m]', 'fontsize', 12);   
axis([0.5 0.7 0 1.25])
plotTickLatex2D

l = legend('J = 21','J = 41','J = 81','J = 161','J = 321','J = 641','J = 1281');


set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
print(gcf,'-depsc2','FMG2.eps')
%% Task 3
exact = dlmread('phi_exact_5000x5000.txt');
%%

figure(2)
clf
set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

plot(0:1/(length(exact)-1):1,exact,0:1/(gridsizes(1)-1):1,phiFMG21(:,halfpoints(1)),0:1/(gridsizes(3)-1):1,phiFMG81(:,halfpoints(3)),0:1/(gridsizes(5)-1):1,phiFMG321(:,halfpoints(5)),0:1/(gridsizes(7)-1):1,phiFMG1281(:,halfpoints(7))) 
y = ylabel('$\Phi$(x,y) [V]', 'fontsize', 12);   
x = xlabel('x [m]', 'fontsize', 12);   
plotTickLatex2D
l = legend('Exact solution','J = 21','J = 81','J = 321','J = 1281');
set(l,'Location','NorthWest')

set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
print(gcf,'-depsc2','FMGN.eps')
clf
plot(0:1/(length(exact)-1):1,exact,0:1/(gridsizes(1)-1):1,phiFMG21(:,halfpoints(1)),0:1/(gridsizes(3)-1):1,phiFMG81(:,halfpoints(3)),0:1/(gridsizes(5)-1):1,phiFMG321(:,halfpoints(5)),0:1/(gridsizes(7)-1):1,phiFMG1281(:,halfpoints(7))) 
y = ylabel('$\Phi$(x,y) [V]', 'fontsize', 12);   
x = xlabel('x [m]', 'fontsize', 12);   
axis([0.5 0.7 0 1.25])
plotTickLatex2D

l = legend('Exact solution','J = 21','J = 81','J = 321','J = 1281');


set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.06, 0]);
print(gcf,'-depsc2','FMGN2.eps')
%% Nice cycle-plots
clc
clf

set(gcf,'renderer','painters','PaperPosition',[0 0 4.7 3]);

V = [81, 41, 21, 11, 21, 41, 81];

W = [81 41 21 11 21 11 21 41 21 11 21 11 21 41 81];

M = [11 21 11 21 41 21 11 21 41 81 41 21 11 21 41 81];


plot(linspace(1,length(V)+1,length(V)), V,'.-','markersize',20);
hold on
plot(linspace(length(V)+3,length(V)+length(W)+3,length(W)), W,'.-','markersize',20);
hold on
plot(linspace(length(W)+3+length(V)+3,length(V)+length(W)+3+3+length(M),length(M)), M,'.-','markersize',20);
hold off

set(gca,'xtick',[])
y = ylabel('Grid size [ ]', 'fontsize', 12);   
%axis([0.5 0.7 0 1.25])
plotTickLatex2D

set(y, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
print(gcf,'-depsc2','cycles.eps')
