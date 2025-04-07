function [Population,FrontNo,CrowdDis] = EnvironmentalSelection_count(Population, N)
    
    [FrontNo,MaxFNo] = countNDSort(Population);

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

function [FrontNo,MaxFNo] = countNDSort(Population)
    % Population: 個体配列 (SOLUTIONクラス配列)
    % 出力: FrontNo(i) = i番目個体の属するFront番号
    PopSize = length(Population);
    CV = ComputeCVM(Population); 
    cv = arrayfun(@(p) sum(max(0,p.cons)), Population);
    Objs = cat(1,Population.objs);

    % "CDPdominates" を用いたペアワイズ比較から被支配数などを求める
    DominateCount = zeros(1,PopSize);
    DominatedSet  = cell(1,PopSize);

    for i = 1 : PopSize
        for j = i+1 : PopSize
            flag = CDPdominates(Objs(i,:), CV(i),cv(i), Objs(j,:), CV(j),cv(j));
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

function flag = CDPdominates(objA, CVA,cvA, objB, CVB,cvB)
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
            if cvA < cvB
                flag = 1;
            elseif cvB < cvA
                flag = -1;
            else
                flag = 0;
            end
        end
    end
end

function CVM = ComputeCVM(Pop)
   % ComputeCVM_Count:
    %   各個体に対し，制約違反となっている制約の“数”をCVMとする
    %
    % 例: CVM(i) = sum_{c=1 to nC} [violation(i,c) > 0 ? 1 : 0]

    cons = cat(1, Pop.cons);
    violation = max(0, cons);

    % violation>0 のboolを合計 => (N×1)
    CVM = sum(violation > 0, 2);
end

