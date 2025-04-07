classdef CM2T_EP < ALGORITHM
% <multi> <real/integer/label/binary/permutation> <constrained>
% Constrained Multi-objective to Two-objective 
% S --- 20 --- The size of subpopulation
% type --- 3 --- The type of aggregation function

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
            [S, type] = Algorithm.ParameterSet(20,3);

            % Parameter 'delta' is the selecting probability from neighbors.
            delta = 0.9;

            %% Generate the weight vectors
            % Parameter 'T' is the size of neighbouring solutions
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            T = ceil(Problem.N/10);

            %% Detect the neighbours of each solution
            B = pdist2(W,W);
            [~,B] = sort(B,2);
            B = B(:,1:T);

            %% Generate random population
            Population = Problem.Initialization();
            feasible_index = find(sum(max(0,Population.cons),2) == 0);
            Nadir = max(Population.objs,[],1);
            Z = min(Population.objs,[],1);

            %% 初期個体を Problem.N 個の初期サブ個体群に収容
            UpdatePopulation = Population;
            Sum_Constraint = sum(max(0,Population.cons),2);
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
                [SubPopulation(i,:),SubPop_Constraint_value(i,:),SubPop_Scalarizing_value(i,:),FrontNo(i,:),CrowdDis(i,:)] = EnvironmentalSelection(UpdatePopulation,Sum_Constraint,Scalarizing_value,S);
            end
        
            Offspring(Problem.N) = SOLUTION;
            Algorithm.EP = Population;
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
                    Nadir = max(Nadir,Offspring(Offspring_gennum).obj);
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
                    [SubPopulation(i,:),SubPop_Constraint_value(i,:),SubPop_Scalarizing_value(i,:),FrontNo(i,:),CrowdDis(i,:)] = EnvironmentalSelection([SubPopulation(i,:),Offspring],Update_Sum_Constraint,Update_Scalarizing_value,S);
                end
                
                %% Output
                
                % サブ個体群を1つの個体群に集約．
                % 出力される個体群サイズは，Problem.NではなくProblem.N*S(重複含む)であることに注意する．
                for i = 1 : Problem.N
                    Population(S*(i-1)+1:S*i) = SubPopulation(i,:);
                end
                % 重複解を削除
                [~,ia,~] = unique(Population.objs,'rows');
                Population = Population(ia);
                Algorithm.EP = updateEP2(Algorithm.EP,Offspring);
            end
        end
    end
end