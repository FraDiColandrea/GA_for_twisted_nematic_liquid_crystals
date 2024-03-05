%% Author Information 
% This function was written by Alicia Sit and Francesco Di Colandrea 
% of the University of Ottawa, 2024.
% See http://arxiv.org/abs/2403.00948 for more details.
%
% LIST OF GA FUNCTIONS
% GA : Main function that will perform the coupled genetic algorithm to 
%      determine integration parameter (beta) and maximum tilt angle 
%      (thetaM) for a twisted nematic liquid crystal cell.
% GAP: Genetic algorithm that solves for thetaM given beta.
% GAB: Genetic algorithm that solves for beta given thetaM.

%% Genetic algorithm beta-thetaM
% INPUT VARIABLES
% om : maximum twist angle of the cell, [0,pi/2].
% V  : externally applied voltage
% iterations: number of generations of the GA
%
% OUTPUT VARIABLES
% bestlistB : list of best beta parameters from each generation
% costlistB : list of cost value for each beta in bestlistB
% bestlistP : list of best thetaM from each generation
% costlistP : list of cost value for each thetaM in bestlistP
function [bestlistB,costlistB,bestlistP,costlistP] = GA(om, V, iterations)
    warning('off','all')
    
    %% Material parameters of LC
    epPar = 12.0;   % dieletric constant parallel to local director
    epPerp = 4.0;   % dieletric constant parallel to local director
    ep0 = 8.85E-12; % vacuum permittivity
    L = 35E-6;      % length of LC cell
    k11 = 6.7E-12;  % splay elastic constant
    k22 = 3.4E-12;  % twist elastic constant
    k33 = 10.6E-12; % bend elastic constant
    
    gamma = (epPar-epPerp)/epPerp;
    alpha = (k33-k22)/k22;
    kappa = (k33-k11)/k11;
    V0 = pi*sqrt(k11/(ep0*(epPar-epPerp)));

    % GA parameters
    probReproduction = 90;
    if V<3 %for low V
        probMutation = 1; 
        generations = 30;
        popSize = 100;  
    else %for high V
        probMutation = 30; 
        generations = 50;
        popSize = 100;  
    end
    
    muMut = 0.0;    % Gaussian noise mean for mutation
    sigmaMut =0.2;  % Gaussion noise standard deviation for mutation 
    elitists = 1;   % # of elistists to carry over from previous generation
    tournamentSelection = popSize - elitists;
    contestants = 4;% # of contestants to compete in tournament
    alph1 = 0.5;    % crossover parameter for parent1
    alph2 = 0.5;    % crossover parameter for parent2

    costlistP = zeros(1,iterations);
    bestlistP = zeros(1,iterations);
    costlistB = zeros(1,iterations);
    bestlistB = zeros(1,iterations);
    
    bet = 0.0;      % initial beta guess

    for n=1:iterations
        [bestlistP(n), costlistP(n)] = GAP(bet);
        phim = bestlistP(n);
        [bestlistB(n), costlistB(n)] = GAB(phim);
        bet = bestlistB(n);
        
        % Plot best value of maximum tilt angle thetaM
        figure(3)
        subplot(1,2,1)
        plot(bestlistP(1:n),'*-')
        title('Best Value ThetaM')
        
        % Plot best value for integration parameter beta
        subplot(1,2,2)
        plot(bestlistB(1:n),'*-')
        title('Best Value Beta')
    
        drawnow
    end



    %% Perform GA for thetaM
    % INPUT VARIABLES:
    % beta : integration parameter
    %
    % OUTPUT VARIABLES:
    % valueBest : best value of thetaM
    % costMin   : minimum cost associated with valueBest
    function [valueBest, costMin] = GAP(beta)
        v = 1;
        pop = randNum(0,pi/2,[popSize,v]); %initial phiM
        elitismParents = zeros(elitists,v);
        bestValue = zeros(generations,v);
        bestCost = zeros(1,generations);
        
        % evalutate cost function
        cost = evalCostP(pop,popSize,beta);
        [minCost,idx] = sort(cost); %get index of lowest cost
        for k=1:elitists
            elitismParents(k,:) = pop(idx(k),:);
        end
        
        % while loop
        j=1;
        while j<=generations
            % do tournament
            winners = zeros(tournamentSelection,v);
            for k=1:tournamentSelection
                winners(k,:) = tournamentP(pop,contestants,beta);
            end
            bestParents = cat(1,elitismParents,winners);
            bestParents = bestParents(randperm(popSize),:); %randomly permute elems
            
            % reproduce
            off1 = zeros(tournamentSelection,v);
            off2 = zeros(tournamentSelection,v);
            for k=1:2:tournamentSelection
                [off1(k,:),off2(k,:)] = reproduceP(bestParents(k,:),bestParents(k+1,:),alph1,alph2,probReproduction);
            end
            offspring = cat(1,off1,off2);
            % mutate
            for k=1:(2*tournamentSelection)
                offspring(k,:) = mutation(offspring(k,:),sigmaMut,muMut,probMutation);
            end
            
            % evalutate cost function
            cost = evalCostP(offspring,popSize,beta);
            [minCost,idx] = sort(cost); %get index of lowest cost
            bestCost(j) = minCost(1);
            bestValue(j,:)= offspring(idx(1),:);
            
            for k=1:length(minCost)
                pop(k,:) = offspring(idx(k),:);
            end
            for k=1:elitists
                elitismParents(k,:) = offspring(idx(k),:);
            end
            %put bestIndividual/bestCost into array
            
            figure(1)
            subplot(1,2,1)
            plot(bestValue(1:j,1),'*-')
            %plot(bestValue(1:j,1),'*-')
            title('Best Value ThetaM')
            
            subplot(1,2,2)
            plot(bestCost(1:j),'*-')
            title('Best Cost ThetaM')
        
            drawnow
        
            j=j+1;
        end
        
        % Best best
        [minCosts,idx] = sort(bestCost); %get index of lowest cost
        bests = zeros(length(minCosts),v);
        for k=1:length(minCosts)
            bests(k,:) = bestValue(idx(k),:);
        end

        valueBest = bests(1,1);
        costMin = minCosts(1);
        %display([num2str(bests(1,1)),', ',num2str(minCosts(1))])
    end
    
    %% Perform GA for beta
    % INPUT VARIABLES:
    % phiM : maximum tilt angle
    %
    % OUTPUT VARIABLES:
    % valueBest : best value of beta
    % costMin   : minimum cost associated with valueBest
    function [valueBest, costMin] = GAB(thetaM)
        v = 1;
        pop = randNum(0,pi/2,[popSize,v]); %initial phiM
        elitismParents = zeros(elitists,v);
        bestValue = zeros(generations,v);
        bestCost = zeros(1,generations);
        
        % evalutate cost function
        cost = evalCostB(pop,popSize,thetaM);
        [minCost,idx] = sort(cost); %get index of lowest cost
        for k=1:elitists
            elitismParents(k,:) = pop(idx(k),:);
        end
        
        % while loop
        j=1;
        while j<=generations
            % do tournament
            winners = zeros(tournamentSelection,v);
            for k=1:tournamentSelection
                winners(k,:) = tournamentB(pop,contestants,thetaM);
            end
            bestParents = cat(1,elitismParents,winners);
            bestParents = bestParents(randperm(popSize),:); %randomly permute elems
            
            % reproduce
            off1 = zeros(tournamentSelection,v);
            off2 = zeros(tournamentSelection,v);
            for k=1:2:tournamentSelection
                [off1(k,:),off2(k,:)] = reproduceB(bestParents(k,:),bestParents(k+1,:),alph1,alph2,probReproduction);
            end
            offspring = cat(1,off1,off2);
            % mutate
            for k=1:(2*tournamentSelection)
                offspring(k,:) = mutation(offspring(k,:),sigmaMut,muMut,probMutation);
            end
            
            % evalutate cost function
            cost = evalCostB(offspring,popSize,thetaM);
            [minCost,idx] = sort(cost); %get index of lowest cost
            bestCost(j) = minCost(1);
            bestValue(j,:)= offspring(idx(1),:);
            
            for k=1:length(minCost)
                pop(k,:) = offspring(idx(k),:);
            end
            for k=1:elitists
                elitismParents(k,:) = offspring(idx(k),:);
            end
            %put bestIndividual/bestCost into array
            
            figure(2)
            subplot(1,2,1)
            plot(bestValue(1:j,1),'*-')
            %plot(bestValue(1:j,1),'*-')
            title('Best Value Beta')
            
            subplot(1,2,2)
            plot(bestCost(1:j),'*-')
            title('Best Cost Beta')
        
            drawnow
        
            j=j+1;
        end
        
        % Best best
        [minCosts,idx] = sort(bestCost); %get index of lowest cost
        bests = zeros(length(minCosts),v);
        for k=1:length(minCosts)
            bests(k,:) = bestValue(idx(k),:);
        end

        valueBest = bests(1,1);
        costMin = minCosts(1);
        %display([num2str(bests(1,1)),', ',num2str(minCosts(1))])
    end
    


    %% Extra Functions

    %%g function
    % theta : tilt angle
    % beta  : integration parameter
    % thetaM: maximum tilt angle
    function g = gFun(theta,beta,thetaM)
        p1 = (sin(thetaM).^2 - sin(theta).^2)/((1+gamma*sin(theta).^2).*(1+gamma*sin(thetaM).^2));
        p2 = (beta.^2*(1+kappa)/(1+alpha)).*(((1+alpha*sin(thetaM).^2).*cos(thetaM).^2).^(-1) - ((1+alpha*sin(theta).^2).*cos(theta).^2).^(-1));
        g = sqrt(p1+p2);
    end
    
    %%Cost function for thetaM
    % pop : population of individuals to determine cost of
    % popSize : size of population
    % beta: integration parameter
    function cost = evalCostP(pop,popSize,beta)
        cost = zeros(1,popSize);
        for i=1:popSize
            phiM = pop(i,:);
            fun = @(x) real(sqrt( 1+kappa.*sin(x).^2 )/( gFun(x,beta,phiM).*(1+gamma.*sin(x).^2 ) ));
            q = integral(fun,0,phiM,'arrayvalued', true);
            cost(i) = ((V/V0)*(pi/2) - q)^2;
        end
    end
    
    %%Cost function for beta
    % pop : population of individuals to determine cost of
    % popSize : size of population
    % thetaM  : maximum tilt angle
    function cost = evalCostB(pop,popSize,thetaM)
        cost = zeros(1,popSize);
        for i=1:popSize
            beta = pop(i,:);
            fun2= @(x) real(beta.*(sqrt(1+kappa.*sin(x).^2)/((gFun(x,beta,thetaM).*cos(x).^2).*(1+alpha.*sin(x).^2))));
            w = integral(fun2,0,thetaM,'arrayvalued', true);
            cost(i) = (0.5*om - w)^2;
        end
    end
    
    % heaviside function
    function h = heaviSide(x)
        if x == 0
            h = 0.5;
        elseif x<0
            h = 0.0;
        elseif x>0
            h = 1.0;
        end
    end
    
    % random number on interval [a,b]
    function r = randNum(a,b,n)
        r = a+(b-a)*rand(n(1),n(2));
    end
    
    %%Reproduction function for thetaM
    % INPUTS
    % p1, p2 : two parents
    % c1, c2 : crossover parameters
    % probReproduction : probability of reproduction
    %
    % OUTPUTS
    % off1, off2 : two offspring from p1, p2
    function [off1, off2] = reproduceP(p1,p2,c1,c2,probReproduction)
        if 100*rand <= probReproduction
            min11 = p1-c1*(p2-p1); min12 = p2-c1*(p1-p2);
            max11 = p2+c1*(p2-p1); max12 = p1+c1*(p1-p2);
            min21 = p1-c2*(p2-p1); min22 = p2-c2*(p1-p2);
            max21 = p2+c2*(p2-p1); max22 = p1+c2*(p1-p2);
            
            norm21 = norm(p2)-norm(p1);
            norm12 = norm(p1)-norm(p2);
    
            off1 = mod(randNum(min11,max11,[1,1])*heaviSide(norm21) + randNum(min12,max12,[1,1])*heaviSide(norm12),pi/2.);
            off2 = mod(randNum(min21,max21,[1,1])*heaviSide(norm21) + randNum(min22,max22,[1,1])*heaviSide(norm12),pi/2);
        else
            off1 = p1;
            off2 = p2;
        end
    end
    
    %%Reproduction function for beta
    % INPUTS
    % p1, p2 : two parents
    % c1, c2 : crossover parameters
    % probReproduction : probability of reproduction
    %
    % OUTPUTS
    % off1, off2 : two offspring from p1, p2
    function [off1, off2] = reproduceB(p1,p2,c1,c2,probReproduction)
        if 100*rand <= probReproduction
            min11 = p1-c1*(p2-p1); min12 = p2-c1*(p1-p2);
            max11 = p2+c1*(p2-p1); max12 = p1+c1*(p1-p2);
            min21 = p1-c2*(p2-p1); min22 = p2-c2*(p1-p2);
            max21 = p2+c2*(p2-p1); max22 = p1+c2*(p1-p2);
            
            norm21 = norm(p2)-norm(p1);
            norm12 = norm(p1)-norm(p2);
    
            off1 = randNum(min11,max11,[1,1])*heaviSide(norm21) + randNum(min12,max12,[1,1])*heaviSide(norm12);
            off2 = randNum(min21,max21,[1,1])*heaviSide(norm21) + randNum(min22,max22,[1,1])*heaviSide(norm12);
        else
            off1 = p1;
            off2 = p2;
        end
    end
    
    %%Mutation function
    % offspring  : the individual to mutate
    % sigma, mu  : Gaussian noise with mean mu and standard deviation sigma
    % proMutation: probability of mutatioin
    function m = mutation(offspring,sigma,mu,probMutation)
        if 100*rand <= probMutation
            m = mod(offspring+(sigma.*randn + mu),pi/2);
            %y = sigma.*randn + mu;
        else
            m = offspring;
        end
    end
    
    %%Tournament for thetaM
    % pop : population of individuals to perform tournament on
    % contestants : # of contestants to compete
    % beta: integration parameter
    function winner = tournamentP(pop,contestants,beta)
        sz =size(pop);
        seed = pop(randperm(sz(1),contestants),:);
        cost = evalCostP(seed,contestants,beta);
    
        [M,I] = min(cost,[],"all"); %get index of lowest cost
        winner = seed(I,:); %winner is the seed with the lowest cost
    end
    
    %%Tournament for beta
    % pop : population of individuals to perform tournament on
    % contestants : # of contestants to compete
    % thetaM: maximum tilt angle
    function winner = tournamentB(pop,contestants,thetaM)
        sz =size(pop);
        seed = pop(randperm(sz(1),contestants),:);
        cost = evalCostB(seed,contestants,thetaM);
    
        [M,I] = min(cost,[],"all"); %get index of lowest cost
        winner = seed(I,:); %winner is the seed with the lowest cost
    end

end