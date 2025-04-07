classdef NSGAIIMo < ALGORITHM
    % <multi> <real/integer/label/binary/permutation> <constrained/none>
    properties
        alpha = 0.2;  % ratio of infeasible solutions
    end

    methods
        function main(Algorithm, Problem)
            %% 1) Generate random population
            pn = Problem.number;
            Population = Problem.Initialization();
            f = 0;
            gen = 1;
            % Initial environment selection
            [Population, FrontNo, CrowdDis] = EnvironmentalSelectionCVM(...
                Population, Problem.N, Algorithm.alpha);
            %filename = sprintf('/Users/azumayuki/Documents/LONs/feasible_ratio_data/RWMOP%d_Mo.csv', pn);
            %filename2 = sprintf('/Users/azumayuki/Documents/LONs/RWMOP%d_Mo_first.csv', pn);

            %f = logFeasibleRatio(Population, gen, filename, f, filename2);
            arch = updateEP2(Population);
            %% Optimization loop
            %% 2) Optimization loop
            while Algorithm.NotTerminated(Population)
                gen = gen + 1;
                % 2.1) Parent selection by Tournament (NSGA-II style)
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);

                % 2.2) Variation (crossover + mutation)
                Offspring  = OperatorGA(Problem, Population(MatingPool));
                arch = updateEP2([arch,Offspring]);
                % 2.3) Environmental selection with alpha ratio
                [Population, FrontNo, CrowdDis] = EnvironmentalSelectionCVM(...
                    [Population, Offspring], Problem.N, Algorithm.alpha);
                %f = logFeasibleRatio(Population, gen, filename, f, filename2);
                if Problem.FE >= Problem.maxFE
                   arch = archive(arch,Problem.N);
                   Population = arch;
                end
            end
        end
    end
end

function [Population, FrontNo, CrowdDis] = EnvironmentalSelectionCVM(PopAll, N, alpha)
    % EnvironmentalSelectionCVM:
    %   1. 全個体に対して CVM を計算し、(f1..fm, CVM) で非優劣ソートを行う。
    %   2. ソートされたリストから実行不可能解を alpha*N 個選択（不足時は全選択）。
    %   3. 残りを実行可能解から選択（不足時は実行不可能解で補完）。
    %   4. 最終集団に対して再度 NDSort と CrowdingDistance を計算する。

    % 1. CVM の計算
    CVM = ComputeCVM(PopAll);

    % 2. (f1..fm, CVM) を用いて非優劣ソート
    PopObjPlus = [cat(1, PopAll.objs), CVM];
    [FrontNoAll, ~] = NDSort(PopObjPlus, length(PopAll));
    CrowdDisAll = CrowdingDistance(PopObjPlus, FrontNoAll);

    % 3. FrontNo と CrowdDis に基づいて全個体をソート (FrontNo 昇順, CrowdDis 降順)
    [~, sortedIdx] = sortrows([FrontNoAll(:), -CrowdDisAll(:)], [1, 2]);
    SortedPop = PopAll(sortedIdx);
    
    % 4. ソートされた個体を実行可能解と実行不可能解に分割
    cons = cat(1, SortedPop.cons);   % (個体数 x 制約数)
    isFeas = all(cons <= 0, 2);      % 全制約を満たすかどうか

    FeasPop = SortedPop(isFeas);     % 実行可能解
    InfePop = SortedPop(~isFeas);    % 実行不可能解

    % 5. 実行不可能解から alpha*N 個選択（不足時は全選択）
    numInfeDesired = round(alpha * N);
    numInfeAvailable = length(InfePop);
    numInfeSel = min(numInfeDesired, numInfeAvailable);

    InfeSel = InfePop(1:numInfeSel);

    % 6. 残りを実行可能解から選択（不足時は実行不可能解で補完）
    numFeasSel = N - numInfeSel;
    numFeasAvailable = length(FeasPop);
    numFeasSel = min(numFeasSel, numFeasAvailable);

    FeasSel = FeasPop(1:numFeasSel);

    % 7. 実行可能解が不足した場合、実行不可能解から補完
    if numFeasSel < (N - numInfeSel)
        shortage = (N - numInfeSel) - numFeasSel;
        if length(InfePop) > numInfeSel
            InfeSel = [InfeSel, InfePop(numInfeSel+1 : min(numInfeSel+shortage, end))];
            FeasSel = FeasSel;  % 変更なし
        end
    end

    % 8. 最終集団の構築
    Population = [InfeSel, FeasSel];

    % 9. 集団が N 未満の場合、ソートされた順に補完
    if length(Population) < N
        remaining = SortedPop(~ismember(SortedPop, Population));
        needed = N - length(Population);
        Population = [Population; remaining(1:min(needed, length(remaining)))];
    end

    % 10. 最終集団が N を超える場合はトリム
    if length(Population) > N
        Population = Population(1:N);
    end

    % 11. 最終集団に対して FrontNo と CrowdDis を再計算
    CVM_final = ComputeCVM(Population);
    PopObjPlusFinal = [cat(1, Population.objs), CVM_final];
    [FrontNo, ~] = NDSort(PopObjPlusFinal, N);
    CrowdDis = CrowdingDistance(PopObjPlusFinal, FrontNo);
end


function CVM = ComputeCVM(Pop)
    % 例: 各制約ごとに違反量を昇順ソートしてRankを付与し、合計をCVMとする
    if isempty(Pop)
        CVM = [];
        return;
    end

    cons = cat(1, Pop.cons);
    [N, nC] = size(cons);

    violation = max(0, cons);  % g_i(x)>0 なら違反
    RankMat = zeros(N, nC);

    for c = 1:nC
        vals = violation(:, c);
        [sortedVal, idx] = sort(vals,'ascend');

        currRank = 0;
        currVal  = -Inf;
        for i = 1:N
            if sortedVal(i) > currVal
                currRank = currRank + 1;
                currVal  = sortedVal(i);
            end
            RankMat(idx(i), c) = currRank - 1; 
            % rank=0 の個体は violation=0 すなわち可行解
        end
    end

    CVM = sum(RankMat, 2);  % (N x 1)
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