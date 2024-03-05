%% Author Information 
% This function was written by Alicia Sit and Francesco Di Colandrea 
% of the University of Ottawa, 2024.
% See http://arxiv.org/abs/2403.00948 for more details.

%% Tilt and Twist GA
% GA to determine the tilt distribution, and then calculate twist 
% distributions.
%
% OUTPUTS
% tiltD  : array, tilt distribution from z/L=[0,1]
% twistD : array, twist distribution from z/L=[0,1]

clear all
close all

% Import list of beta and thetaM found from GA.m
betaThetaM = importdata([pwd,'\betaThetaM.txt']);

phiM = pi/4; % max twist angle (rad)

%%
i=1; % Choose voltage according to betaThetaM index

beta = betaThetaM(i,2);
thetaM = betaThetaM(i,3);
delta = 0.07;

tic
[tiltList,tiltCost] = GAtilt(delta,beta,thetaM);
twistList = twist(beta,thetaM,tiltList);

tiltD = cat(2,tiltList,thetaM,flip(tiltList));
twistD = cat(2,twistList,0.5*phiM,phiM-flip(twistList));
toc

%% Plot tilt and twist distributions
figure(3)
z = (0:1/20:1);
yyaxis left
plot(z,tiltD*180/pi,'r-*')
axis([0 1 0 90])
ylabel(['Tilt angle (',char(176),')'])
xlabel('z/L')
title(['Applied voltage V = ',num2str(betaThetaM(i,1)),' V'])

yyaxis right
plot(z,twistD*180/pi,'b-*')
axis([0 1 0 phiM*180/pi])
ylabel(['Twist angle (',char(176),')'])
legend('Tilt distribution, \theta(z)','Twist distribution, \phi(z)','Location','southeast')
set(get(gca,'ylabel'),'rotation',-90,'VerticalAlignment','bottom')

ax = gca;
ax.YAxis(1).Color = 'r';
ax.YAxis(2).Color = 'b';

