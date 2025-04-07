classdef NSGAIIeps < ALGORITHM
    % <multi> <real/integer/label/binary/permutation> <constrained/none>
    % NSGA-II with Epsilon-Constraint Handling
    %   初期εは「初期個体群内の 0.05N番目に大きいCV」を使用
    methods
        function main(Algorithm,Problem)
            pn = Problem.number;
            %% -- 1. 初期集団の生成 --
            Population = Problem.Initialization();  % サイズN
            
            % 集団サイズ
            N = length(Population);

            % ---- 1.1 各個体の制約違反量CVを計算 ----
            CV = arrayfun(@(p) sum(max(0,p.cons)), Population);

            % ---- 1.2 CVを大きい順にソートし、0.05N番目の値を初期εとする ----
            sortedCV_desc = sort(CV,'descend'); 
            idx = ceil(0.05*N); % 切り上げ or floor() で調整可能
            if idx < 1
                idx = 1; % 念のため
            elseif idx > N
                idx = N;
            end
            epsilon_0 = sortedCV_desc(idx);

            % ---- 1.3 環境選択を実行 (Generation=0時点) ----
            f = 0;
            gen = 1;
            epsilon    = epsilon_0;        % 今のε
            Tc               = 0.9*ceil(Problem.maxFE/Problem.N);
            cp               = 2;
            [Population,FrontNo,CrowdDis] = EnvironmentalSelectionEpsilon(Population, Problem.N, epsilon);

            filename = sprintf('/Users/azumayuki/Documents/LONs/feasible_ratio_data/RWMOP%d_eps.csv', pn);
            filename2 = sprintf('/Users/azumayuki/Documents/LONs/RWMOP%d_eps_first.csv', pn);

            f = logFeasibleRatio(Population, gen, filename, f, filename2);
            arch = updateEP2(Population);
            %% -- 2. 世代ループ --
            while Algorithm.NotTerminated(Population)
                gen = gen + 1;

                % 2.1 現在の世代に応じて ε を更新 (線形例)
                

                % 2.2 トーナメント選択
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);

                % 2.3 交叉・突然変異
                Offspring  = OperatorGA(Problem, Population(MatingPool));
                arch = updateEP2([arch,Offspring]);
                % 2.4 環境選択（ε制約法により可行解優先）
                [Population,FrontNo,CrowdDis] = EnvironmentalSelectionEpsilon([Population,Offspring], Problem.N, epsilon);
                if gen < Tc
                    epsilon = epsilon_0 * ((1 - (gen / Tc)) ^ cp);
                else
                    epsilon = 0;
                end
                % 2.5 可行解割合をログ
                f = logFeasibleRatio(Population, gen, filename, f, filename2);
                if Problem.FE >= Problem.maxFE
                   arch = archive(arch,Problem.N*2);
                   Population = arch;
                end
            end
        end
    end
end

%%========= 下記サブルーチンたち ==========%%
function [Population,FrontNo,CrowdDis] = EnvironmentalSelectionEpsilon(Population, N, epsilon)
   
    CV = arrayfun(@(p) sum(max(0,p.cons)), Population);
    FeasibleMask = (CV <= epsilon);
    FeasiblePop  = Population(FeasibleMask);
    InfeasiblePop= Population(~FeasibleMask);
    CV_infe = CV(~FeasibleMask);

    % 可行解 -> NDsort + CrowdingDistance
    if ~isempty(FeasiblePop)
        [FrontNoF,MaxFNoF] = NDSort(cat(1,FeasiblePop.obj),length(FeasiblePop));
        CrowdDisF = CrowdingDistance(cat(1,FeasiblePop.obj),FrontNoF);
    else
        FrontNoF = [];
        CrowdDisF= [];
        MaxFNoF  = 0;
    end

    % 不可行解 -> まとめて1つのフロントにしてCV小さい順でソート(例)
    if ~isempty(InfeasiblePop)
        [sortedCV,idx_infe] = sort(CV_infe,'ascend');
        InfeasiblePop = InfeasiblePop(idx_infe);
        FrontNoI  = MaxFNoF + 1 * ones(1,length(InfeasiblePop));
        CrowdDisI = zeros(1,length(InfeasiblePop));
    else
        FrontNoI = [];
        CrowdDisI= [];
    end

    % 結合
    MergedPop   = [FeasiblePop, InfeasiblePop];
    MergedFNo   = [FrontNoF,   FrontNoI];
    MergedCDis  = [CrowdDisF,  CrowdDisI];
    MaxFNo      = max([0, MergedFNo]);

    % NSGA-II と同様、Frontが小さい順に埋めてN個を選ぶ
    Next = false(1,length(MergedPop));
    for f = 1:MaxFNo
        idxF = find(MergedFNo==f);
        if sum(Next) + length(idxF) <= N
            Next(idxF) = true;
        else
            [~,rankC] = sort(MergedCDis(idxF),'descend');
            remain = N - sum(Next);
            Next(idxF(rankC(1:remain))) = true;
            break;
        end
    end

    Population = MergedPop(Next);
    FrontNo    = MergedFNo(Next);
    CrowdDis   = MergedCDis(Next);
end

function t = logFeasibleRatio(Pop, gen, filename,f ,filename2)
    CV = arrayfun(@(p) sum(max(0,p.cons)), Pop);
    feasibleCount = sum(CV == 0);
    
    ratio = feasibleCount / length(Pop);
    if f == 0
        if ratio > 0
            f = 1;
            fid = fopen(filename2,'a');
            fprintf(fid,"%d,%f\n",gen);
            fclose(fid);
        end
    end
    t = f;
    fid = fopen(filename,'a');
    fprintf(fid,"%d,%f\n",gen,ratio);
    fclose(fid);
end