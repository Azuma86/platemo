classdef PM_SR < ALGORITHM
    % <multi> <real/integer/label/binary/permutation> <constrained/none>
    methods
        function main(Algorithm,Problem)
            pn = Problem.number;
            rng(125)
            N = Problem.N;
            D = Problem.D;
            
            lhsSample = lhsdesign(N, D);  % Generate LHS sample
            lowerBounds = repmat(Problem.lower, N, 1);
            upperBounds = repmat(Problem.upper, N, 1);
            PopDec = lowerBounds + lhsSample .* (upperBounds - lowerBounds);
            Population = Problem.Evaluation(PopDec);
            proM = 0.9;
            disM = 20;
            child = 10;
            feaGen = zeros(1,N);
            stableCount = zeros(1,N);
            BestCon = zeros(1,N);
            stableThreshold = 20;
            for i = 1:N
                BestCon(i) = sum(max(0,Population(i).cons),2);
                if BestCon(i)  == 0
                    stableCount(i) = stableThreshold;
                end   
            end
            % CSVファイルの初期化: 決定変数と目的関数、制約値をカラム分割
            csvFilename = sprintf('/Users/azumayuki/Documents/LONs/data09-20-pre/RWMOP%d_SR.csv',pn);
            %csvFilename = sprintf(csvFilename,proM,disM);
            objCount = length(Population(1).objs);
            conCount = length(Population(1).cons);
            decHeader = arrayfun(@(x) sprintf('X_%d',x), 1:D, 'UniformOutput', false);
            objHeader = arrayfun(@(x) sprintf('Obj_%d',x), 1:objCount, 'UniformOutput', false);
            conHeader = arrayfun(@(x) sprintf('Con_%d',x), 1:conCount, 'UniformOutput', false);
            header = [{'Gen', 'ID'}, decHeader{:}, objHeader{:}, conHeader{:}];
            
            fid = fopen(csvFilename, 'w');
            fprintf(fid, '%s,', header{1:end-1});
            fprintf(fid, '%s\n', header{end});
            fclose(fid);
            currentGeneration = 0;
            while Algorithm.NotTerminated(Population)
                currentGeneration = currentGeneration + 1;

                % 現世代個体の情報をCSVへ
                fid = fopen(csvFilename, 'a');
                for i = 1:length(Population)
                    individual = Population(i);
                    dec = individual.decs;
                    x = [];
                    for t = 1:length(x)
                        dec(x(t)) = round(dec(x(t)));
                    end
                    decStr = strjoin(arrayfun(@num2str, dec, 'UniformOutput', false), ',');
                    objStr = strjoin(arrayfun(@num2str, individual.objs, 'UniformOutput', false), ',');
                    conStr = strjoin(arrayfun(@num2str, individual.cons, 'UniformOutput', false), ',');
                    fprintf(fid, '%d,%d,%s,%s,%s\n', currentGeneration, i, decStr, objStr, conStr);
                end
                fclose(fid);

                % 次世代生成
                Offspring = Population;  % Offspringは後で上書きする
                for i = 1:N
                    Parent = Population(i);
                    if stableCount(i) < stableThreshold
                        Candidates = repmat(Parent,1,child+1);
                        for j = 1:child
                            mutatedDec = polynomialMutation(Parent.decs, Problem.lower, Problem.upper, proM, disM);
                            Candidates(j+1) = Problem.Evaluation(mutatedDec);
                        end
                        Offspring(i) = SelectByStochasticRankingWithRankCV(Candidates, 0.45);
                        if BestCon(i) <= sum(max(0,Offspring(i).cons),2)%一定世代間制約違反量が改善されなければ終了
                            stableCount(i) = stableCount(i) + 1;
                        else
                            stableCount(i) = 0;
                            BestCon(i) = sum(max(0,Offspring(i).cons),2);
                        end
                        if sum(max(0,Offspring(i).cons),2)==0%実行可能解になれば終了
                            feaGen(i) = currentGeneration;
                            stableCount(i) = stableThreshold;
                        end
                    else
                        % 安定状態なら突然変異なしでそのまま
                        Offspring(i) = Parent;
                    end
                end

                % 選択は上で既に11個体中から最良を選んでいるので、そのまま
                Population = Offspring;
       
                if all(stableCount >= stableThreshold)
                    CV = arrayfun(@(p) sum(max(0,p.cons)), Population);
                    feasibleCount = sum(CV == 0);
                    ratio = feasibleCount / N;
                    mean_gen = mean(feaGen(CV==0));
                    fprintf("実行可能解割合：%f 平均到達世代:%f\n",ratio,mean_gen)
                    break;
                end
            end
        end
    end
end

function mutatedOffspring = polynomialMutation(Offspring, lower, upper, proM, disM)
    D = length(Offspring);
    for i = 1:D
        if rand < proM
            mu = rand;
            yl = lower(i);
            yu = upper(i);
            if mu <= 0.5
                deltaq = ((2*mu) + (1 - 2*mu)*((1-(Offspring(i)-yl)/(yu-yl))^(disM+1)))^(1/(disM+1)) - 1;
                Offspring(i) = Offspring(i) + deltaq*(yu-yl);
            else
                deltaq = 1 - ((2*(1-mu)) + (2*(mu-0.5)*((1-(yu-Offspring(i))/(yu-yl))^(disM+1))))^(1/(disM+1));
                Offspring(i) = Offspring(i) + deltaq*(yu-yl);
            end
        end
    end
    mutatedOffspring = min(max(Offspring, lower), upper);
end

function bestInd = SelectByStochasticRankingWithRankCV(popArray, pf)
    % popArray: 親 + 子 個体配列(INDIVIDUAL)
    % pf: SRで"目的関数比較"を行う確率 (ex: 0.45)
    %
    %  1) まず全個体のFrontNo(ランク)をNDSortで求める
    %  2) CV(制約違反量)を計算
    %  3) SRバブルソート: 目的比較の際は (ランク <=> CV)で優劣を決める
    
    N = length(popArray);
    
    % 1. 非優劣ソート (すべての目的が最小化想定)
    PopObj = cat(1, popArray.objs);
    FrontNo = NDSort(PopObj, N);  % => 1,2,3..., ∞

    % 2. 制約違反量
    CV = zeros(N,1);
    for i = 1:N
        CV(i) = sum(max(0, popArray(i).cons));
    end

    % 3. SRバブルソート
    maxPass = 2*N;
    pass = 1;
    swapped = true;
    while pass <= maxPass && swapped
        swapped = false;
        for i = 1:(N-1)
            if ~compareSR(popArray(i), popArray(i+1), ...
                          FrontNo(i), FrontNo(i+1), ...
                          CV(i), CV(i+1), pf)
                % swap
                [popArray(i), popArray(i+1)] = deal(popArray(i+1), popArray(i));
                [FrontNo(i), FrontNo(i+1)] = deal(FrontNo(i+1), FrontNo(i));
                [CV(i), CV(i+1)]           = deal(CV(i+1), CV(i));
                swapped = true;
            end
        end
        pass = pass + 1;
    end
    
    bestInd = popArray(1);  % ソート後の先頭個体を返す
end

function better = compareSR(indA, indB, rankA, rankB, cvA, cvB, pf)
    % SRによる隣接要素比較
    % 戻り値 better = true => "AはBより優れている"
    %          false => "Bの方が優れている"
    %
    % Step:
    %   1) 確率pfで "目的関数比較(=ランク+CV)"
    %   2) 確率(1-pf)で "制約違反量"
    
    mu = rand();
    if mu <= pf
        %% 目的関数比較 => まずランク，ランク同じならCV
        if rankA < rankB
            better = true;
        elseif rankA > rankB
            better = false;
        else
   %{
            %同ランク => CV小さい方を上位
            if cvA < cvB
                better = true;
            elseif cvA > cvB
                better = false;
            else
                % CVも同じ => 適当にAを上位とする(同等)
   %}
                better = true;
            %end
        end
    else
        %% 制約比較 => CVが小さい方が優先
        better = (cvA < cvB);
    end
end