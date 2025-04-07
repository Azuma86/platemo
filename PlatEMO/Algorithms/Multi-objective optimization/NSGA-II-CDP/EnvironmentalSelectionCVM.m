function [Population,FrontNo,CrowdDis] = EnvironmentalSelectionCVM(Population, N)
    
    [FrontNo,MaxFNo] = CVMNDSort(Population);

    CrowdDis = CrowdingDistance(Population.objs,FrontNo);

    Next = false(1,length(Population));

    for f = 1 : MaxFNo
        inds = find(FrontNo==f);
        if length(inds)+sum(Next) <= N
            Next(inds) = true;
        else
            [~,Rank]   = sort(CrowdDis(inds),'descend');
            Next(inds(Rank(1:(N-sum(Next))))) = true;
            break;
        end
    end

    % 4. 出力
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end

function [FrontNo,MaxFNo] = CVMNDSort(Population)
    % Population: 個体配列 (SOLUTIONクラス配列)
    % 出力: FrontNo(i) = i番目個体の属するFront番号
    PopSize = length(Population);
    CV = ComputeCVM(Population);  % 制約違反量 (0が可行)
    Objs = cat(1,Population.objs);

    % "CDPdominates" を用いたペアワイズ比較から被支配数などを求める
    DominateCount = zeros(1,PopSize);
    DominatedSet  = cell(1,PopSize);

    for i = 1 : PopSize
        for j = i+1 : PopSize
            flag = CDPdominates(Objs(i,:), CV(i), Objs(j,:), CV(j));
            if flag == 1
                % iがjを支配
                DominatedSet{i} = [DominatedSet{i}, j];
                DominateCount(j) = DominateCount(j) + 1;
            elseif flag == -1
                % jがiを支配
                DominatedSet{j} = [DominatedSet{j}, i];
                DominateCount(i) = DominateCount(i) + 1;
            end
        end
    end

    % 各Frontに振り分け
    FrontNo = inf(1,PopSize);
    front = 1;
    current = find(DominateCount==0);
    while ~isempty(current)
        FrontNo(current) = front;
        nextTier = [];
        for i = current
            for j = DominatedSet{i}
                DominateCount(j) = DominateCount(j) - 1;
                if DominateCount(j) == 0
                    nextTier = [nextTier, j];
                end
            end
        end
        front = front + 1;
        current = nextTier;
    end
    MaxFNo = front - 1;
end

function flag = CDPdominates(objA, CVA, objB, CVB)
    % 1) Aが可行でB不可行 => A支配
    % 2) 両者可行 => 通常の多目的支配
    % 3) 両者不可行 => 制約違反が小さい方が支配
    if CVA == 0 && CVB > 0
        flag = 1;
    elseif CVA > 0 && CVB == 0
        flag = -1;
    elseif CVA == 0 && CVB == 0
        % 両者可行 -> 多目的支配
        if all(objA <= objB) && any(objA < objB)
            flag = 1;
        elseif all(objB <= objA) && any(objB < objA)
            flag = -1;
        else
            flag = 0;
        end
    else
        % 両者不可行 -> CVが小さい方が優位
        if CVA < CVB
            flag = 1;
        elseif CVB < CVA
            flag = -1;
        else
            flag = 0;
        end
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


