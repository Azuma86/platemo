classdef PM_Mo < ALGORITHM
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
            feaGen = zeros(1,N);
            stableCount = zeros(1,N);
            BestCon = zeros(1,N);
            stableThreshold = 20;
            child = 10;
            for i = 1:N
                BestCon(i) = sum(max(0,Population(i).cons),2);
                if BestCon(i)  == 0
                    stableCount(i) = stableThreshold;
                end   
            end

            % CSVファイルの初期化: 決定変数と目的関数、制約値をカラム分割
            csvFilename = sprintf('/Users/azumayuki/Documents/LONs/data09-20-pre/RWMOP%d_Mo.csv',pn);
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
                        % 親個体から10個の子個体を生成
                        Candidates = repmat(Parent,1,child+1);
                        for j = 1:child
                            mutatedDec = polynomialMutation(Parent.decs, Problem.lower, Problem.upper, proM, disM);
                            Candidates(j+1) = Problem.Evaluation(mutatedDec);
                        end
                        Offspring(i) = SelectByObjCVM(Candidates);
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

function bestIndividual = SelectByObjCVM(popArray)
    % popArray: 親+子 (INDIVIDUAL配列, 例えばサイズ=11)
    % 1) Compute CVM for each individual
    % 2) NDSort( [obj, CVM] ), then crowding distance
    % 3) 「可行解があれば最良可行解を選択」, なければ順位1位を選択

    N = length(popArray);
    % (1) 目的関数値
    PopObj = cat(1, popArray.objs);

    % (2) CVM 計算 (制約違反量をrank化して合計)
    CVM = ComputeCVM(popArray);

    % (3) [PopObj, CVM] で非優劣ソート
    CombinedObj = [PopObj, CVM];
    FrontNo = NDSort(CombinedObj, N);
    idxSet = find(FrontNo == 1);
    popArray = popArray(idxSet);
    if isscalar(idxSet)
        bestIndividual = popArray;
    else
        for k = 1 : length(idxSet)
            % 総制約違反量を計算
            cvVal = sum(max(0, popArray(k).cons));
            if cvVal == 0
                % 可行解が見つかったら即座にそれを返す
                bestIndividual = popArray(k);
                return;
            end
        end
        r = randi(length(idxSet));      
        bestIndividual = popArray(r);
    end
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

