classdef PPS_CM2T < ALGORITHM
% <multi/many> <real/integer> <constrained>
% Push and pull search algorithm
% S --- 20 --- The size of subpopulation
% nr    ---   2 --- Maximum number of solutions replaced by each offspring
% type --- 1 --- The type of aggregation function
% rate_update_weight --- 0.05 --- Ratio of updated weight vectors
% wag                ---  250 --- Iteration interval of utilizing AWA
%------------------------------- Reference --------------------------------
% Z. Fan, W. Li, X. Cai, H. Li, C. Wei, Q. Zhang, K. Deb, and E. Goodman,
% Push and pull search for solving constrained multi-objective optimization
% problems, Swarm and Evolutionary Computation, 2019, 44(2): 665-679.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Wenji Li

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [S,nr,type,rate_update_weight,wag] = Algorithm.ParameterSet(20,2,1,0.05,250);

            %% Generate the weight vectors
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            % Transformation on W
            W = 1./W./repmat(sum(1./W,2),1,size(W,2));
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            Z          = min(Population.objs,[],1);

            %% Evaluate the Population
            Tc               = 0.9 * ceil(Problem.maxFE/Problem.N);
            last_gen         = 20;
            search_stage     = 1; % 1 for push stage,otherwise,it is in pull stage.
            change_threshold = 0;
            epsilon_k        = 0;
            epsilon_0        = 0;
            cp               = 2;
            alpha            = 0.95;
            tao              = 0.05;
            max_change       = 1;
            ideal_points     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            nadir_points     = zeros(ceil(Problem.maxFE/Problem.N),Problem.M);
            update = 0;
            arch             = archive(Population,Problem.N);
            EP = Population;
            % Size of external elite
            nEP = ceil(Problem.N*1.5);
            nus = rate_update_weight*Problem.N;
            Pi         = ones(Problem.N,1);
            oldObj     = max(abs((Population.objs-repmat(Z,Problem.N,1)).*W),[],2);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                gen        = ceil(Problem.FE/Problem.N);
                pop_cons   = Population.cons;
                cv         = overall_cv(pop_cons);
                population = [Population.decs,Population.objs,cv];
                rf         = sum(cv <= 1e-6) / Problem.N;
                ideal_points(gen,:) = Z;
                nadir_points(gen,:) = max(population(:,Problem.D + 1 : Problem.D + Problem.M),[],1);
                
                % The maximumrate of change of ideal and nadir points rk is calculated.
                if gen >= last_gen
                    max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen);
                end
                % The value of e(k) and the search strategy are set.
                if gen < Tc
                    if Problem.FE > Problem.maxFE*0.5 && search_stage == 1
                        search_stage = -1;
                        epsilon_0 = max(population(:,end),[],1);
                        epsilon_k = epsilon_0;
                    end
                    if max_change <= change_threshold && search_stage == 1
                        search_stage = -1;
                        epsilon_0 = max(population(:,end),[],1);
                        epsilon_k = epsilon_0;
                    end
                   
                    if search_stage == -1
                        epsilon_k = update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,gen,Tc,cp);
                    end
                else
                    epsilon_k = 0;
                end
                
                if search_stage == 1 % Push Stage
                    if ~mod(ceil(Problem.FE/Problem.N),10)
                        % Allocation of computing resources
                        newObj    = max(abs((Population.objs-repmat(Z,Problem.N,1)).*W),[],2);
                        DELTA     = (oldObj-newObj)./oldObj;
                        Temp      = DELTA <= 0.001;
                        Pi(~Temp) = 1;
                        Pi(Temp)  = (0.95+0.05*DELTA(Temp)/0.001).*Pi(Temp);
                        oldObj    = newObj;
                    end
                    for subgeneration = 1 : 5
                        % Choose I
                        Bounday = find(sum(W<1e-3,2)==1)';
                        I = [Bounday,TournamentSelection(10,floor(Problem.N/5)-length(Bounday),-Pi)];
    
                        % Evolve each solution in I
                        Offspring(1:length(I)) = SOLUTION();
                        for i = 1 : length(I)
                            % Choose the parents
                            if rand < 0.9
                                P = B(I(i),randperm(size(B,2)));
                            else
                                P = randperm(Problem.N);
                            end
    
                            % Generate an offspring
                            Offspring(i) = OperatorGAhalf(Problem,Population(P(1:2)));
    
                            % Update the ideal point
                            Z = min(Z,Offspring(i).obj);

                            switch type
                                case 1
                                    % Tchebycheff approach
                                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
                                    g_new = max(repmat(abs(Offspring(i).obj-Z),length(P),1).*W(P,:),[],2);
                                case 2
                                    % Tchebycheff approach with normalization
                                    Zmax  = max(Population.objs,[],1);
                                    g_old = max(abs(Population(P).objs-repmat(Z,length(P),1))./repmat(Zmax-Z,length(P),1).*W(P,:),[],2);
                                    g_new = max(repmat(abs(Offspring(i).obj-Z)./(Zmax-Z),length(P),1).*W(P,:),[],2);
                            end
                            Population(P(find(g_old>=g_new,nr))) = Offspring(i);
                        end
                    end
                else  % Pull Stage  &&  An improved epsilon constraint-handling is employed to deal with constraints
                    if update == 0
                        UpdatePopulation = Population;
                        Sum_Constraint = sum(max(0,Population.cons),2);
                        for i = 1 : Problem.N
                            % Update the ideal point
                            switch type
                                case 1
                                    % Tchebycheff approach
                                    Zmax  = max(Population.objs,[],1);
                                    Scalarizing_value = max(abs(Population.objs-repmat(Z,Problem.N,1)).*W(i,:),[],2);
                                case 2
                                    % Tchebycheff approach with normalization
                                    Zmax  = max(Population.objs,[],1);
                                    Scalarizing_value = max(abs(Population.objs-repmat(Z,Problem.N,1))./repmat(Zmax - Z,Problem.N,1).*W(i,:),[],2);
                            end
                        
                            % スカラー化関数と制約違反量を目的関数値とする個体群を作成．
                            [SubPopulation(i,:),SubPop_Constraint_value(i,:),SubPop_Scalarizing_value(i,:)] = EpsilonEnvironmentalSelection(UpdatePopulation,Sum_Constraint,Scalarizing_value,S,epsilon_k);
                        end
                        update = update + 1;
                        Offspring(Problem.N) = SOLUTION;
                    end
                    Offspring_gennum = 1;
                    [Minimize_Constraint_value,~] = min(SubPop_Constraint_value, [], 2);
                    % Compute r
                    opt_index = ones(Problem.N,1);
                    r = ones(Problem.N,1);
                    for i = 1 : Problem.N
                        minimize_constraint_index_insubpop = find(SubPop_Constraint_value(i,:) == Minimize_Constraint_value(i));
                        if length(minimize_constraint_index_insubpop) > 1
                            [~,index] = min(SubPop_Scalarizing_value(i,minimize_constraint_index_insubpop));
                            opt_index(i) = index(1);
                        else
                            opt_index(i) = minimize_constraint_index_insubpop;
                        end
                        r(i) = length( find(SubPop_Scalarizing_value(i,:)<=SubPop_Scalarizing_value(i,opt_index(i))) ) / S;
                    end
                    I = randperm(Problem.N);
                    for i = 1 : Problem.N
                        % 注目するサブ個体群No
                        SubPopNo = I(i);
                        
                        % Choose the parent subpopulation
                        if rand < delta
                            P = B(I(i),randperm(end));
                        else
                            P = randperm(Problem.N);
                        end
                        
                        % Choose the parent solutions
                        ParentSolutionDecs(3) = SOLUTION;
                        if rand < r(SubPopNo)
                            ParentSolutionDecs(1) = SubPopulation(P(1),opt_index(P(1)));
                        else
                            ParentSolutionDecs(1) = SubPopulation(P(1),randsample(S,1));
                        end
                        decs = SubPopulation(B(I(i),:),:);
                        decs = decs(:);
                        random_order = randperm(numel(decs(:)));
                        ParentSolutionDecs(2) = decs(random_order(1));
                        ParentSolutionDecs(3) = decs(random_order(2));
                        
                        % Generate offspring using DE
                        Offspring(Offspring_gennum) = OperatorDE(Problem,ParentSolutionDecs(1),ParentSolutionDecs(2),ParentSolutionDecs(3));
                        
                        % Update reference point
                        Zmax = max(Zmax,Offspring(Offspring_gennum).obj);
                        Z = min(Z,Offspring(Offspring_gennum).obj);
                        Offspring_gennum = Offspring_gennum + 1;
                    end
                    
                    % 現在の現個体群と全ての子個体群を合わせた Gloal.N + S 個でサブ個体群の更新を行う．
                    Offspring_Cons = Offspring.cons;
                    Offspring_Objs = Offspring.objs;
                    
                    for i = 1 : Problem.N
                        % 子個体に目的関数（スカラー化関数と制約違反量）をセット．
                        Offspring_Sum_Constraint = sum(max(0,Offspring_Cons),2);
                        Curent_Sum_Constraint = sum(max(0,SubPopulation(i,:).cons),2);
                        Curent_Objs = SubPopulation(i,:).objs;
                        switch type
                            case 1
                                % Tchebycheff approach
                                Offspring_Scalarizing_value = max(abs(Offspring_Objs-repmat(Z,Problem.N,1)).*W(i,:),[],2);
                                Curent_Scalarizing_value = max(abs(Curent_Objs-repmat(Z,S,1)).*W(i,:),[],2);
                            case 2
                                % Tchebycheff approach with normalization
                                Offspring_Scalarizing_value = max(abs(Offspring_Objs-repmat(Z,Problem.N,1))./repmat(Zmax-Z,Problem.N,1).*W(i,:),[],2);
                                Curent_Scalarizing_value = max(abs(Curent_Objs-repmat(Z,S,1))./repmat(Zmax-Z,S,1).*W(i,:),[],2);    
                        end
                        Update_Sum_Constraint = [Curent_Sum_Constraint; Offspring_Sum_Constraint];
                        Update_Scalarizing_value = [Curent_Scalarizing_value; Offspring_Scalarizing_value];
                        
                        % Gloal.N + S 個でサブ個体群 i の更新．
                        [SubPopulation(i,:),SubPop_Constraint_value(i,:),SubPop_Scalarizing_value(i,:)] = EpsilonEnvironmentalSelection([SubPopulation(i,:),Offspring],Update_Sum_Constraint,Update_Scalarizing_value,S,epsilon_k);
                    end
          
                    for i = 1 : Problem.N
                        Population(S*(i-1)+1:S*i) = SubPopulation(i,:);
                    end

                    if ~mod(ceil(Problem.FE/Problem.N),wag/5) && ~isempty(EP)
                        [~,ia,~] = unique(EP.objs,'rows');
                        EP = EP(ia);
                        [Population,Del,W] = updateWeightAWA2(Population,SubPopulation,W,Z,EP,nus,T);
                        n = length(find(Del));
                        %% 初期個体を Problem.N 個の初期サブ個体群に収容
                        UpdatePopulation = Population;
                        Sum_Constraint = sum(max(0,Population.cons),2);
                        for i = 1 : n
                            switch type
                                case 1
                                    % Tchebycheff approach
                                    Zmax  = max(Population.objs,[],1);
                                    Scalarizing_value = max(abs(Population.objs-repmat(Z,Problem.N,1)).*W(i,:),[],2);
                                case 2
                                    % Tchebycheff approach with normalization
                                    Zmax  = max(Population.objs,[],1);
                                    Scalarizing_value = max(abs(Population.objs-repmat(Z,Problem.N,1))./repmat(Zmax - Z,Problem.N,1).*W(i,:),[],2);                            
                            end
                            % スカラー化関数と制約違反量を目的関数値とする個体群を作成．
                            [SubPopulation(Problem.N-n+i,:),SubPop_Constraint_value(Problem.N-n+i,:),SubPop_Scalarizing_value(Problem.N-n+i,:)] = EpsilonEnvironmentalSelection(UpdatePopulation,Sum_Constraint,Scalarizing_value,S,epsilon_k);
                        end
                    end
              
                    % サブ個体群を1つの個体群に集約．
                    % 出力される個体群サイズは，Problem.NではなくProblem.N*S(重複含む)であることに注意する．
                    for i = 1 : Problem.N
                        Population(S*(i-1)+1:S*i) = SubPopulation(i,:);
                    end
                    % 重複解を削除
                    [~,ia,~] = unique(Population.objs,'rows');
                    Population = Population(ia);
                   
                end
                EP = updateEPAWA(EP,Offspring,nEP);
                
                % Output the non-dominated and feasible solutions.
                arch = archive([arch,Population],Problem.N);
                if Problem.FE >= Problem.maxFE
                    Population = arch;
                end
                
            end
        end
    end
end

% Calculate the Maximum Rate of Change
function max_change = calc_maxchange(ideal_points,nadir_points,gen,last_gen)
    delta_value = 1e-6 * ones(1,size(ideal_points,2));
    rz = abs((ideal_points(gen,:) - ideal_points(gen - last_gen + 1,:)) ./ max(ideal_points(gen - last_gen + 1,:),delta_value));
    nrz = abs((nadir_points(gen,:) - nadir_points(gen - last_gen + 1,:)) ./ max(nadir_points(gen - last_gen + 1,:),delta_value));
    max_change = max([rz, nrz]);
end

% The Overall Constraint Violation
function result = overall_cv(cv)
    cv(cv <= 0) = 0;cv = abs(cv);
    result = sum(cv,2);
end

function result = update_epsilon(tao,epsilon_k,epsilon_0,rf,alpha,gen,Tc,cp)
    if rf < alpha
        result = (1 - tao) * epsilon_k;
    else
        result = epsilon_0 * ((1 - (gen / Tc)) ^ cp);
    end
end