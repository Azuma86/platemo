classdef CM2TAWARS < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Constrained Multi-objective to Two-objective 
% S --- 20 --- The size of subpopulation
% type --- 3 --- The type of aggregation function
% rate_update_weight --- 0.05 --- Ratio of updated weight vectors
% wag                ---  100 --- Iteration interval of utilizing AWA
% rate_evol          ---  0.8 --- Ratio of iterations to evolve with only MOEA/D
% rate_evol2 ---  0.5 --- Ratio of iterations to evolve
%evaluation --- 1 ---final population of EP



%------------------------------- Reference --------------------------------
% Takafumi Fukase, Naoki Masuyama, Yuusuke Nojima, Hisao Ishibuchi, 
% A Constrained Multi objective Evolutionary Algorithm Based on 
% Transformation to Two objective Optimization Problems, In Proc. of
% FAN2019, Toyama, 2019.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
           %% Parameter setting
            % Parameter 'type' is the type of scalarizing function.
            % Parameter 'S' is the size of subpopulation.
            [S, type,rate_update_weight,wag, rate_evol,rate_evol2,evaluation] = Algorithm.ParameterSet(20,3,0.05,100,0.8,0.5,1);

            % Parameter 'delta' is the selecting probability from neighbors.
            delta = 0.9;
            % number of update vector
            nus = rate_update_weight*Problem.N;
            %% Generate the weight vectors
            % Parameter 'T' is the size of neighbouring solutions
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);
             % Size of external elite
            nEP = ceil(Problem.N*1.5);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            feasible_index = find(sum(max(0,Population.cons),2) == 0);
            Nadir = max(Population.objs,[],1);
            Z = min(Population.objs,[],1);
            epsilon = Inf;
            EP = Population;
            %% 初期個体を Problem.N 個の初期サブ個体群に収容
            UpdatePopulation = Population;
            Sum_Constraint = sum(max(0,Population.cons),2);
            arch = updateEP2(Population);
            for i = 1 : Problem.N
                switch type
                    case 1
                        % PBI approach
                        normW   = sqrt(sum(W(i,:).^2,2));
                        normP   = sqrt(sum((Population.objs-repmat(Z,Problem.N,1)).^2,2));
                        CosineP = sum((Population.objs-repmat(Z,Problem.N,1)).*W(i,:),2)./normW./normP;
                        Scalarizing_value   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                    case 2
                        % Tchebycheff approach with normalization
                        Scalarizing_value = max(abs(Population.objs-repmat(Z,Problem.N,1))./repmat(Nadir - Z,Problem.N,1).*W(i,:),[],2);
                    case 3
                        % Tchebycheff approach
                        Scalarizing_value = max(abs(Population.objs-repmat(Z,Problem.N,1)).*W(i,:),[],2);
                end
                % スカラー化関数と制約違反量を目的関数値とする個体群を作成．
                [SubPopulation(i,:),SubPop_Constraint_value(i,:),SubPop_Scalarizing_value(i,:)] = EpsilonEnvironmentalSelection(UpdatePopulation,Sum_Constraint,Scalarizing_value,S,epsilon);
            end
        
            Offspring(Problem.N) = SOLUTION;
           
            %% Optimization
            while Algorithm.NotTerminated(Population)
                % サブ個体群内で子個体の生成と参照点の更新．
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
                for i = 1 : Problem.N
                    % 注目するサブ個体群No                    
                    % Choose the parent subpopulation
                    if rand < delta
                        P = B(i,randperm(end));
                    else
                        P = randperm(Problem.N);
                    end
                    
                    % Choose the parent solutions
                    ParentSolutionDecs(3) = SOLUTION;
                    if rand < r
                        ParentSolutionDecs(1) = SubPopulation(P(1),opt_index(P(1)));
                    else
                        ParentSolutionDecs(1) = SubPopulation(P(1),randsample(S,1));
                    end
                    decs = SubPopulation(B(i,:),:);
                    decs = decs(:);
                    random_order = randperm(numel(decs(:)));
                    ParentSolutionDecs(2) = decs(random_order(1));
                    ParentSolutionDecs(3) = decs(random_order(2));
                    
                    % Generate offspring using DE
                    Offspring(Offspring_gennum) = OperatorDE(Problem,ParentSolutionDecs(1),ParentSolutionDecs(2),ParentSolutionDecs(3));
                    
                    % Update reference point
                    Nadir = max(Nadir,Offspring(Offspring_gennum).obj);
                    Z = min(Z,Offspring(Offspring_gennum).obj);
                    Offspring_gennum = Offspring_gennum + 1;
                end
                
                % 現在の現個体群と全ての子個体群を合わせた Gloal.N + S 個でサブ個体群の更新を行う．
                Offspring_Cons = Offspring.cons;
                Offspring_Objs = Offspring.objs;
                arch = updateEP2([arch,Offspring]);
                for i = 1 : Problem.N
                    % 子個体に目的関数（スカラー化関数と制約違反量）をセット．
                    Offspring_Sum_Constraint = sum(max(0,Offspring_Cons),2);
                    Curent_Sum_Constraint = sum(max(0,SubPopulation(i,:).cons),2);
                    Curent_Objs = SubPopulation(i,:).objs;
                    switch type
                        case 1
                            % PBI approach
                            normW   = sqrt(sum(W(i,:).^2,2));
                            normP   = sqrt(sum((Offspring_Objs-repmat(Z,Problem.N,1)).^2,2));
                            normP_Curent   = sqrt(sum((Curent_Objs-repmat(Z,S,1)).^2,2));
                            CosineP = sum((Offspring_Objs-repmat(Z,Problem.N,1)).*W(i,:),2)./normW./normP;
                            CosineP_Curent = sum((Curent_Objs-repmat(Z,S,1)).*W(i,:),2)./normW./normP_Curent;
                            Offspring_Scalarizing_value   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                            Curent_Scalarizing_value   = normP_Curent.*CosineP_Curent + 5*normP_Curent.*sqrt(1-CosineP_Curent.^2);
                        case 2
                            % Tchebycheff approach with normalization
                            Offspring_Scalarizing_value = max(abs(Offspring_Objs-repmat(Z,Problem.N,1))./repmat(Nadir-Z,Problem.N,1).*W(i,:),[],2);
                            Curent_Scalarizing_value = max(abs(Curent_Objs-repmat(Z,S,1))./repmat(Nadir-Z,S,1).*W(i,:),[],2);
                        case 3
                            % Tchebycheff approach
                            Offspring_Scalarizing_value = max(abs(Offspring_Objs-repmat(Z,Problem.N,1)).*W(i,:),[],2);
                            Curent_Scalarizing_value = max(abs(Curent_Objs-repmat(Z,S,1)).*W(i,:),[],2);
                    end
                    
                    Update_Sum_Constraint = [Curent_Sum_Constraint; Offspring_Sum_Constraint];
                    Update_Scalarizing_value = [Curent_Scalarizing_value; Offspring_Scalarizing_value];
                    
                    % Gloal.N + S 個でサブ個体群 i の更新．
                    [SubPopulation(i,:),SubPop_Constraint_value(i,:),SubPop_Scalarizing_value(i,:)] = EpsilonEnvironmentalSelection([SubPopulation(i,:),Offspring],Update_Sum_Constraint,Update_Scalarizing_value,S,epsilon);
                end
                
                %% Output
                                %% Output
                EP = updateEPAWA(EP,Offspring,nEP);
                if Problem.FE > rate_evol*Problem.maxFE
                    % サブ個体群を1つの個体群に集約．
                    % 出力される個体群サイズは，Problem.NではなくProblem.N*S(重複含む)であることに注意する．
                    for i = 1 : Problem.N
                        Population(S*(i-1)+1:S*i) = SubPopulation(i,:);
                    end
                    if ~mod(ceil(Problem.FE/Problem.N),wag/5) && ~isempty(EP)
                        [~,ia,~] = unique(EP.objs,'rows');
                        EP = EP(ia);
                        [Population,SubPopulation,SubPop_Constraint_value,SubPop_Scalarizing_value,Del,W] = updateWeightAWA3(Population,SubPopulation,SubPop_Constraint_value,SubPop_Scalarizing_value,W,Z,EP,nus,T);
                        n = length(find(Del));
                        %% 初期個体を Problem.N 個の初期サブ個体群に収容
                        UpdatePopulation = Population;
                        Sum_Constraint = sum(max(0,Population.cons),2);
                        for i = 1 : n
                            switch type
                                case 1
                                    % PBI approach
                                    normW   = sqrt(sum(W(i,:).^2,2));
                                    normP   = sqrt(sum((Population.objs-repmat(Z,Problem.N,1)).^2,2));
                                    CosineP = sum((Population.objs-repmat(Z,Problem.N,1)).*W(i,:),2)./normW./normP;
                                    Scalarizing_value   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                                case 2
                                    % Tchebycheff approach with normalization
                                    Scalarizing_value = max(abs(Population.objs-repmat(Z,length(Population),1))./repmat(Nadir - Z,length(Population),1).*W(Problem.N-n+i,:),[],2);
                                case 3
                                    % Tchebycheff approach
                                    Scalarizing_value = max(abs(Population.objs-repmat(Z,length(Population),1)).*W(Problem.N-n+i,:),[],2);
                            end
                            % スカラー化関数と制約違反量を目的関数値とする個体群を作成．
                            [SubPopulation(Problem.N-n+i,:),SubPop_Constraint_value(Problem.N-n+i,:),SubPop_Scalarizing_value(Problem.N-n+i,:)] = EpsilonEnvironmentalSelection(UpdatePopulation,Sum_Constraint,Scalarizing_value,S,epsilon);
                        end
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
                
                %epsilon 更新
                if Problem.FE > rate_evol2*Problem.maxFE
                    epsilon = maxCV*((Problem.maxFE-Problem.FE)/(Problem.maxFE*(1-rate_evol2)));
                else
                    maxCV = max(sum(max(0,Population.cons),2));
                end
                % Output the non-dominated and feasible solutions.
                if Problem.FE >= Problem.maxFE
                    if evaluation == 2
                        Population = arch;
                    end
                end
            end
        end
    end
end