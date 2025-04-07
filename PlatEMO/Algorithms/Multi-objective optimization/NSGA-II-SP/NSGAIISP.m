classdef NSGAIISP < ALGORITHM
    % <multi> <real/integer/label/binary/permutation> <constrained/none>
    % NSGA-II with SP-based (Tessema & Yen) constraint handling, 
    % assuming that the problem class outputs .cons for all constraints.
    %
    % Reference:
    %  - K. Deb, et al., "A fast and elitist multiobjective genetic algorithm: NSGA-II."
    %  - T. Tessema and G. Yen, "A Self Adaptive Penalty Function Based on 
    %    Stochastic Ranking for Constrained Evolutionary Optimization."
    methods
        function main(Algorithm, Problem)
            pn = Problem.number;
            Population = Problem.Initialization(); 
            f = 0;
            [Population, FrontNo, CrowdDis] = EnvironmentalSelectionSP(Population, Problem.N);
            gen = 1;
            filename = sprintf('/Users/azumayuki/Documents/LONs/feasible_ratio_data/RWMOP%d_SP.csv', pn);
            filename2 = sprintf('/Users/azumayuki/Documents/LONs/RWMOP%d_SP_first.csv', pn);

            f = logFeasibleRatio(Population, gen, filename, f, filename2);
            arch = updateEP2(Population);
            while Algorithm.NotTerminated(Population)
                gen = gen+1;
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);
                
                Offspring = OperatorGA(Problem, Population(MatingPool));
                arch = updateEP2([arch,Offspring]);
                disp(length(arch))
                [Population, FrontNo, CrowdDis] = EnvironmentalSelectionSP([Population,Offspring], Problem.N);
                f = logFeasibleRatio(Population, gen, filename, f, filename2);
                if Problem.FE >= Problem.maxFE
                   arch = archive(arch,Problem.N*2);
                   Population = arch;
                end
            end
        end
    end
end


function [Population, FrontNo, CrowdDis] = EnvironmentalSelectionSP(Population, N)
% 環境選択：
% 1) 個体ごとの修正目的値 (SP) を計算
% 2) 修正目的値に基づき非優越ソート
% 3) 混み合い距離 (CrowdingDistance) を計算
% 4) 上位N個体を選択

    % ---------- 1. SPによる修正目的値を計算 ----------
    PopObjModified = CalcModifiedObjectivesSP(Population);

    % ---------- 2. 非優越ソート (修正目的値で) ----------
    % 下記NDsortはPlatEMO標準関数を想定
    [FrontNo, MaxFNo] = NDSort(PopObjModified, N); 
    Next = FrontNo < MaxFNo; % 完全に選ばれたフロント
    
    % ---------- 3. 混み合い距離の計算 (修正目的値で) ----------
    CrowdDis = CrowdingDistance(PopObjModified, FrontNo);

    % ---------- 4. 最終フロントの補充 ----------
    Last = find(FrontNo == MaxFNo);
    [~,Rank] = sort(CrowdDis(Last), 'descend');
    Selected = Rank(1 : N - sum(Next));

    Next(Last(Selected)) = true;
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end

function PopObjModified = CalcModifiedObjectivesSP(Population)
% SP(距離測定値+適応型ペナルティ)を用いて
% 修正目的値 \tilde{f}_i'(x) を計算し、行列として返す
%
% 入力:
%  - Population : 個体配列 (SOLUTIONクラスの配列)
% 出力:
%  - PopObjModified (個体数 x M) : 各個体・各目的の修正目的値

    PopObj = cat(1, Population.objs);  % (個体数 x M)
    PopCon = cat(1, Population.cons); 
    PopCon(PopCon < 0) = 0;
    [N, M] = size(PopObj);

    % (1) 各目的の [min, max] を取得 (正規化のため)
    f_min = min(PopObj, [], 1);
    f_max = max(PopObj, [], 1);
    c_max = max(PopCon, [], 1);
    c_max(c_max==0) = 1e-12; 
    denom = f_max - f_min;
    denom(denom==0) = 1e-12;  % 定数目的を避けるための小さい値
    
    PopObjNorm = (PopObj - f_min) ./ denom;  % (N x M)
    PopConNorm = PopCon ./ c_max;  % (N x C)
    v = mean(PopConNorm, 2);        % (N x 1)
    v1 = repmat(v,1,M);
    numFeasible = sum(v == 0);  
    eta = numFeasible / N;
  
   
    % (4) 距離測定値 D_i(x):
    %     実行可能(CV=0)なら正規化目的値, 非実行可能(CV>0)ならCV
    D = zeros(N, M);
    for i = 1 : N
        if eta == 0
            D(i,:) = v1(i,:);
        else
            D(i,:) = sqrt( (PopObjNorm(i,:)).^2 + (v1(i,:)).^2 );
        end
    end

    % (5) ペナルティの計算 P_{1,i}(x), P_{2,i}(x)
    P1 = zeros(N, M);
    P2 = zeros(N, M);
    for i = 1 : N
        if eta == 0
            P1(i,:) = 0;  
        else
            P1(i,:) = v1(i,:);
        end
        if v(i)==0
            P2(i,:) = 0;
        else
            P2(i,:) = PopObjNorm(i,:);
        end

    end

    % (7) SPによる修正目的値 \tilde{f}_i'(x) = D_i(x) + η P_{1,i}(x) + (1-η) P_{2,i}(x)
    PopObjModified = D + eta .* P2 + (1 - eta) .* P1;  % (N x M)
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