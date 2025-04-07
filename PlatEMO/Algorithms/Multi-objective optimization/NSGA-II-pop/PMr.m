classdef PMr < ALGORITHM
    % <multi> <real/integer/label/binary/permutation> <constrained/none>
    % Nondominated sorting genetic algorithm II
    
    %------------------------------- Reference --------------------------------
    % K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
    % multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
    % Evolutionary Computation, 2002, 6(2): 182-197.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    methods
        function main(Algorithm,Problem)
            rng(125)
            N = Problem.N;
            D = Problem.D;
            lhsSample = lhsdesign(N, D);  % Generate LHS sample
            lowerBounds = repmat(Problem.lower, N, 1);
            upperBounds = repmat(Problem.upper, N, 1);
            PopDec = lowerBounds + lhsSample .* (upperBounds - lowerBounds);
            Population = Problem.Evaluation(PopDec);
            proM = 1.0;
            disM = 20;
            t = 1;
            %% Generate random population
            %Population = Problem.Initialization();
            % CSVファイルを開き、ヘッダーを書き込む
            csvFilename = '/Users/azumayuki/Documents/LONs/p_data.csv';
            header = {'run', 'fit1', 'node1','x1_1','x2_1', 'fit2','node2','x1_1','x2_1'};
            %header = {'Gen', 'ID','Con'};
            fid = fopen(csvFilename, 'w');
            fprintf(fid, '%s %s %s %s %s %s %s %s %s\n', header{:});
            %fprintf(fid, '%s,%s,%s\n', header{:});
            fclose(fid);
            
            while Algorithm.NotTerminated(Population)
                % 現時点の世代（Generation）情報を取得、または管理
                gen = ceil(Problem.FE/Problem.N);
                p = Population;
                % polyminal mutaion
                % Offspringの初期化
                Offspring = Population;
                % 各親個体に多項式変異を適用
                for i = 1:Problem.N
                    Parent = Population(i);
                    lowerBound = Problem.lower;
                    upperBound = Problem.upper;
                    
                    %polynomialMutation関数から変異した子を生成
                    Offspring(i) = Problem.Evaluation(polynomialMutation(Parent, lowerBound, upperBound, proM, disM));
                    CV_pop = sum(max(0,Population(i).cons),2);
                    CV_off = sum(max(0,Offspring(i).cons),2);

                    if CV_pop > CV_off
                        Population(i) = Offspring(i);
                    end
                end
                % 各個体の遺伝子と目的関数の値を保存
                fid = fopen(csvFilename, 'a'); % ファイルを追記モードで開く
                for i = 1:N
                    CV1 = sum(max(0,p(i).cons),2);
                    CV2 = sum(max(0,Population(i).cons),2);
                    x1_1 = p(i).dec(1);
                    x2_1 = p(i).dec(2);
                    x1_2 = Population(i).dec(1);
                    x2_2 = Population(i).dec(2);
                    t1 = sprintf('%d_%d', gen-1,i);
                    t2 = sprintf('%d_%d', gen,i);
                    % データをCSVに書き込む
                    fprintf(fid,'%d %d %s %d %d %d %s %d %d\n',t,CV1,t1,x1_1,x2_1,CV2,t2,x1_2,x2_2);
                end
                fclose(fid);
            end
        end
    end
end


function mutatedOffspring = polynomialMutation(parent, lower, upper, proM, disM)
    % Extract decision variables from the parent
    Offspring = parent.decs;
    D = length(Offspring);
    mutatedValues = Offspring;
    for i = 1:D
        if rand < proM / D
            mu = rand;
            yl = lower(i);
            yu = upper(i);
            x = Offspring(i); % Use the decision variable
            if mu <= 0.5
                Offspring       = min(max(Offspring,lower),upper);
                Offspring(i) = Offspring(i)+(yu-yl).*((2.*mu+(1-2.*mu).*...
                                  (1-(Offspring(i)-yl)./(yu-yl)).^(disM+1)).^(1/(disM+1))-1);    
            else
                Offspring(i) = Offspring(i)+(yu-yl).*(1-(2.*(1-mu)+2.*(mu-0.5).*...
                      (1-(yu-Offspring(i))./(yu-yl)).^(disM+1)).^(1/(disM+1)));
            end
        end
    end
    mutatedOffspring = Offspring;
end