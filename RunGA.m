%% Author Information 
% This function was written by Alicia Sit and Francesco Di Colandrea 
% of the University of Ottawa, 2024.
% See http://arxiv.org/abs/2403.00948 for more details.

%% Running GA.m
% Main file to execute the coupled GA to find the integration parameter 
% (beta) and maximum tilt angle (thetaM) for a given maximum twist angle 
% (om) and externally applied voltage (V)
%
% OUPUTS
% beta  : array of numerically calculated beta for a given om, at each V
% thetaM: array of numerically calculated thetaM for a given om, at each V

clear all
close all

fileName = [pwd,'\betaThetaM.txt']; % file to save calculations to

iter = 10;  % # of generations to run GA
om = pi/4;  % set maximum twist angle
Vmin = 6.0; % minimum voltage
Vmax = 7.0; % maximum voltage
Vstep = 0.5;% voltage increment
V = (Vmin:Vstep:Vmax)./(2.0*sqrt(2.0));
%V = 6.0;

beta = zeros(length(V),1);     
thetaM = zeros(length(V),1);   
costBeta = zeros(length(V),1);
costPhiM = zeros(length(V),1);

betaThetaM = zeros(length(V),3);
betaThetaM(:,1) = V';

for i=1:length(V)
    display(num2str(length(V)-i))
    
    % GA for best beta-thetaM
    tic
    [bestlistB,costlistB,bestlistP,costlistP] = GA(om,V(i),iter);
    toc
    
    % best values and costs
    [~, idx] = sort(costlistB);
    if idx(1)==1
        minCostB = costlistB(idx(2));
        minCostP = costlistP(idx(2));
    
        bestB = bestlistB(idx(2));
        bestP = bestlistP(idx(2));
    
%         display([num2str(bestB),', ',num2str(minCostB)])
%         display([num2str(bestP),', ',num2str(minCostP)])
    else
        minCostB = costlistB(idx(1));
        minCostP = costlistP(idx(1));
    
        bestB = bestlistB(idx(1));
        bestP = bestlistP(idx(1));
    
%         display([num2str(bestB),', ',num2str(minCostB)])
%         display([num2str(bestP),', ',num2str(minCostP)])

    end
    
    beta(i) = bestB;
    thetaM(i) = bestP;
    costBeta(i) = minCostB;
    costPhiM(i) = minCostP;

end

betaThetaM(:,2) = beta';
betaThetaM(:,3) = thetaM';

writematrix(betaThetaM,fileName,'Delimiter','tab')

