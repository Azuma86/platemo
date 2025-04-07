function [FrontNo,MaxFNo] = CVNDSort(Population)
    % Population: 個体配列 (SOLUTIONクラス配列)
    % 出力: FrontNo(i) = i番目個体の属するFront番号
    PopSize = length(Population);
    CV = arrayfun(@(p) sum(max(0,p.cons)), Population);  % 制約違反量 (0が可行)
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