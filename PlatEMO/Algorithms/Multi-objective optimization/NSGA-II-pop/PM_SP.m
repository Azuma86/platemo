classdef PM_SP < ALGORITHM
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
            csvFilename = sprintf('/Users/azumayuki/Documents/LONs/data09-20-pre/RWMOP%d_SP.csv',pn);
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
                    %{
                    dec(5) = round(dec(5));
                    dec(6) = round(dec(6));
                    dec(4) = round(dec(4));
                    dec(7) = round(dec(7));
                    %}
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
                        
                        Offspring(i) = SelectBySP(Candidates);
                        if BestCon(i) <= sum(max(0,Offspring(i).cons),2)%一定世代間制約違反量が改善されなければ終了
                            stableCount(i) = stableCount(i) + 1;
                        else
                            stableCount(i) = 0;
                            BestCon(i) = sum(max(0,Offspring(i).cons),2);
                        end
                        if sum(max(0,Offspring(i).cons),2)==0%実行可能解になれば終了
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

function bestInd = SelectBySP(popArray)
    % popArray : 1×(child+1) の個体配列 (親+子)
    % Step:
    %   1) 修正目的値 PopObjModified を計算 (N×M) 
    %   2) 非優劣ソート + crowdingDistance
    %   3) rank1, crowdDist大から順に1個取り => bestInd
    %
    % 参考: EnvironmentalSelectionSPを簡略化

    N = length(popArray);
    PopObj = cat(1, popArray.objs);   % (N×M)
    PopCon = cat(1, popArray.cons);   % (N×C)

    % -- SPで "修正目的値" を計算 (N×M)
    PopObjModified = CalcModifiedObjectivesSP(PopObj, PopCon);

    % -- 非優劣ソート
    FrontNo = NDSort(PopObjModified, N);
    idxSet = find(FrontNo == 1);
    popArray = popArray(idxSet);
    if isscalar(idxSet)
        bestInd = popArray;
    else
        for k = 1 : length(idxSet)
            % 総制約違反量を計算
            cvVal = sum(max(0, popArray(k).cons));
            if cvVal == 0
                % 可行解が見つかったら即座にそれを返す
                bestInd = popArray(k);
                return;
            end
        end
        r = randi(length(idxSet));      
        bestInd = popArray(r);
    end
end

function PopObjModified = CalcModifiedObjectivesSP(PopObj, PopCon)
    % PopObj: (N×M)  目的関数
    % PopCon: (N×C)  制約 g_i(x) <= 0 (正の部分だけ違反量)
    % Tessema & Yenの式に準拠 (例)
    
    [N, M] = size(PopObj);
    PopCon(PopCon<0) = 0;  % 不要部分を0
    violation = PopCon;    % max(0,PopCon)後

    % (1) 目的値正規化
    f_min = min(PopObj,[],1);
    f_max = max(PopObj,[],1);
    denom_f = f_max - f_min;
    denom_f(denom_f==0) = 1e-12;
    PopObjNorm = (PopObj - f_min)./denom_f;

    % (2) 制約違反正規化
    c_max = max(violation, [], 1);
    c_max(c_max==0) = 1e-12;
    PopConNorm = violation ./ c_max;  % (N×C)

    v = mean(PopConNorm,2);  % 各個体の平均違反
    numFeasible = sum(v==0);
    eta = 0;   % 可行率

    % (3) 距離測定 D
    D = zeros(N,M);
    for i=1:N
        if eta==0
            % 全不可行なら => v(i)を使う (ここではベクトルでなく一様に v(i)
            D(i,:) = v(i);
        else
            % sqrt( obj^2 + v^2 )
            D(i,:) = sqrt(PopObjNorm(i,:).^2 + v(i)^2);
        end
    end

    % (4) Penalty P1, P2
    P1 = zeros(N,M);
    P2 = zeros(N,M);
    for i=1:N
        if eta>0
            P1(i,:) = v(i);
        end
        if v(i)>0
            P2(i,:) = PopObjNorm(i,:);
        end
    end

    % (5) \tilde{f}'(x) = D + eta * P2 + (1-eta)*P1
    PopObjModified = D + eta.*P2 + (1-eta).*P1;
end