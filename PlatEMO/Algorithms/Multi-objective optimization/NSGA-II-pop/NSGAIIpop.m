classdef NSGAIIpop < ALGORITHM
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
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
            % CSVファイルを開き、ヘッダーを書き込む
            csvFilename = '/Users/azumayuki/Documents/LONs/population_data.csv';
            header = {'Gen', 'ID', 'X', 'Obj','Con'};
            %header = {'Gen', 'ID','Con'};
            fid = fopen(csvFilename, 'w');
            fprintf(fid, '%s,%s,%s,%s,%s\n', header{:});
            %fprintf(fid, '%s,%s,%s\n', header{:});
            fclose(fid);
            
           while Algorithm.NotTerminated(Population)
                % 現時点の世代（Generation）情報を取得、または管理
                currentGeneration = ceil(Problem.FE/Problem.N);
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));

                % 各個体の遺伝子と目的関数の値を保存
                fid = fopen(csvFilename, 'a'); % ファイルを追記モードで開く
                for i = 1:length(Population)
                    individual = Population(i);
                    genes = strjoin(arrayfun(@num2str, individual.decs, 'UniformOutput', false), ','); % 遺伝子情報を文字列に変換
                    objectives = strjoin(arrayfun(@num2str, individual.objs, 'UniformOutput', false), ','); % 目的関数値を文字列に変換
                    CV = sum(max(0,Population(i).cons),2);
                    % データをCSVに書き込む
                    fprintf(fid, '%d,%d,"%s","%s","%d"\n', currentGeneration, i, genes, objectives,CV);
                    %fprintf(fid, '%d,%d,"%d"\n', currentGeneration, i, genes, objectives,CV);
                end
                fclose(fid);
            
                % 次の個体群を生成
                [Population, FrontNo, CrowdDis] = EnvironmentalSelection([Population, Offspring], Problem.N);
            end
        end
    end
end