classdef NSGAIIchev < ALGORITHM
    % <multi> <real/integer/label/binary/permutation> <constrained/none>
    % NSGA-II with Constrained Domination Principle (CDP)
    %
    % This implementation stores the ratio of feasible solutions each generation. 
    methods
        function main(Algorithm, Problem)
            pn = Problem.number;
            %% 1. 初期集団の生成
            Population = Problem.Initialization();  % ランダム初期化
            % CDPによる環境選択でFrontNo,CrowdDisを得る
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection_chev(Population, Problem.N);
            f = 0;  
            gen = 1;
            %% 1.1 CSVファイルを用意（世代ごとの可行解割合を記録）
            % （必要に応じてファイル名や出力形式を変更してください）
            filename = sprintf('/Users/azumayuki/Documents/LONs/feasible_ratio_CV/RWMOP%d_chev.csv', pn);
            filename2 = sprintf('/Users/azumayuki/Documents/LONs/feasible_first_CV/RWMOP%d_chev_first.csv', pn);
            %f = logFeasibleRatio(Population, gen, filename, f, filename2);
            arch = updateEP2(Population);
            %% 2. メインループ
            while Algorithm.NotTerminated(Population)
                gen = gen + 1;

                % 2.1 トーナメント選択 (CDP考慮)
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);

                % 2.2 交叉・突然変異
                Offspring  = OperatorGA(Problem, Population(MatingPool));
                arch = updateEP2([arch,Offspring]);
                % 2.3 次世代選択（CDPベース）
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection_chev([Population,Offspring], Problem.N);

                %f = logFeasibleRatio(Population, gen, filename, f, filename2);
                if Problem.FE >= Problem.maxFE
                   arch = archive(arch,Problem.N*2);
                   Population = arch;
                end
            end
        end
    end
end

function t = logFeasibleRatio(Pop, gen, filename,f ,filename2)
    CV = sum(max(0,Pop.cons),2);
    disp(length(find(CV==0))/length(Pop))
    ratio = length(find(CV==0)) / length(Pop);
    if f == 0
        if ratio > 0
            f = 1;
            fid = fopen(filename2,'a');
            fprintf(fid,"%d,%f\n",gen);
            fclose(fid);
        end
    end
    t = f;
    disp(ratio)

    fid = fopen(filename,'a');
    fprintf(fid,"%d,%f\n",gen,ratio);
    fclose(fid);

end