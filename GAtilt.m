%% Author Information 
% This function was written by Alicia Sit and Francesco Di Colandrea 
% of the University of Ottawa, 2024.
% See http://arxiv.org/abs/2403.00948 for more details.

%% Genetic algorithm to solve the tilt distribution theta(z)
% INPUT VARIABLES
% delta : uniform noise to perturb solution
% beta  : integration parameter
% thetaM  : maximum tilt angle
%
% OUTPUT VARIABLES
% bestlistB : list of best beta parameters from each generation
% costlistB : list of cost value for each beta in bestlistB
% bestlistP : list of best thetaM from each generation
% costlistP : list of cost value for each thetaM in bestlistP
function [tiltList, tiltCost] = GAtilt(delta,beta,thetaM)
    warning('off','all')
    
    %% Material parameters of LC
    epPar = 12.0;   % dieletric constant parallel to local director
    epPerp = 4.0;   % dieletric constant parallel to local director
    ep0 = 8.85E-12; % vacuum permittivity
    k11 = 6.7E-12;  % splay elastic constant
    k22 = 3.4E-12;  % twist elastic constant
    k33 = 10.6E-12; % bend elastic constant
    
    gamma = (epPar-epPerp)/epPerp;
    alpha = (k33-k22)/k22;
    kappa = (k33-k11)/k11;
    
    L = 35E-6;              % length of LC cell
    z = (0:L/20:L/2)./L;
    % the algorithm can't compute at L/2, so shorten span
    z = z(1:(length(z)-1)); 
    
    % GA parameters
    probReproduction = 90;
    probMutation = 1;
    generations = 20;
    popSize = 50;  
    
    muMut = 0.0;
    sigmaMut =0.2;
    elitists = 1;
    tournamentSelection = popSize - elitists;
    contestants = 4;
    alph1 = 0.5;
    alph2 = 0.5;
    
    fun2 = @(x) real(sqrt(1+kappa.*sin(x).^2)./gFun(x,beta,thetaM));
    normInt = integral(fun2,0,thetaM,'arrayvalued', true);
    
    tiltList = zeros(1,length(z));
    tiltCost = zeros(1,length(z));
    prevTilt = 0;
    
    for r=1:length(z)
        [tiltList(r),tiltCost(r)] = GA(z(r),r,prevTilt);
        prevTilt = tiltList(r); 
        
        figure(2)
        plot(z,tiltList*180/pi,'*-')
        title('Tilt Distribution, \theta(z)')
        xlabel('z/L')
        ylabel(['Tilt angle (',char(176),')'])
        drawnow
    end
       
    
    %% GA
    % INPUTS
    % d : distance
    % iter : current iteration
    % prevTilt : tilt angle from previous iteration
    function [valueBest,costMin] = GA(d,iter,prevTilt)
        v = 1;
        if iter==0
            pop = randNum(0,pi/2,[popSize,v]); %initial phiZ
        else%if iter<90
            pop = prevTilt+randNum(0,delta,[popSize,v]);
    %     elseif iter>=90
    %         pop = prevTilt+randNum(0,0.01,[popSize,v]);
        end
        elitismParents = zeros(elitists,v);
        bestValue = zeros(generations,v);
        bestCost = zeros(1,generations);
    
        % evalutate cost function
        cost = evalCost(pop,popSize,d,beta,thetaM);
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
                winners(k,:) = tournament(pop,contestants,d,beta,thetaM);
            end
            bestParents = cat(1,elitismParents,winners);
            bestParents = bestParents(randperm(popSize),:); %randomly permute elems
    
            % reproduce
            off1 = zeros(tournamentSelection,v);
            off2 = zeros(tournamentSelection,v);
            for k=1:2:tournamentSelection
                [off1(k,:),off2(k,:)] = reproduce(bestParents(k,:),bestParents(k+1,:),alph1,alph2,probReproduction);
            end
            offspring = cat(1,off1,off2);
            % mutate
            for k=1:(2*tournamentSelection)
                offspring(k,:) = mutation(offspring(k,:),sigmaMut,muMut,probMutation);
            end
    
            % evalutate cost function
            cost = evalCost(offspring,popSize,d,beta,thetaM);
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
            plot(bestValue(1:j,1)*180/pi,'*-')
            %plot(bestValue(1:j,1),'*-')
            title('Best Value ThetaM')
            xlabel('Generation')
    
            subplot(1,2,2)
            plot(bestCost(1:j),'*-')
            title('Best Cost ThetaM')
            xlabel('Generation')
    
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
    
    
    %% Functions
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
    
    % g function
    function g = gFun(theta,beta,thetaM)
        p1 = (sin(thetaM).^2 - sin(theta).^2)/((1+gamma*sin(theta).^2).*(1+gamma*sin(thetaM).^2));
        p2 = (beta.^2*(1+kappa)/(1+alpha)).*(((1+alpha*sin(thetaM).^2).*cos(thetaM).^2).^(-1) - ((1+alpha*sin(theta).^2).*cos(theta).^2).^(-1));
        g = sqrt(p1+p2);
    end
    
    % cost function for phiM
    function cost = evalCost(pop,popSize,d,beta,thetaM)
        cost = zeros(1,popSize);
        for i=1:popSize
            tilt = pop(i,:);
            fun = @(x) real(sqrt(1+kappa.*sin(x).^2)./gFun(x,beta,thetaM));
            q = integral(fun,0,tilt,'arrayvalued', true);
            cost(i) = (2*d*normInt - q)^2;
        end
    end
    
    % reproduction function for theta
    function [off1, off2] = reproduce(p1,p2,c1,c2,probReproduction)
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
    
    % mutation function
    function m = mutation(offspring,sigma,mu,probMutation)
        if 100*rand <= probMutation
            m = mod(offspring+(sigma.*randn + mu),pi/2);
            %y = sigma.*randn + mu;
        else
            m = offspring;
        end
    end
    
    % tournament for theta
    function winner = tournament(pop,contestants,d,beta,thetaM)
        sz =size(pop);
        seed = pop(randperm(sz(1),contestants),:);
        cost = evalCost(seed,contestants,d,beta,thetaM);
    
        [M,I] = min(cost); %get index of lowest cost
        winner = seed(I,:); %winner is the seed with the lowest cost
    end

end