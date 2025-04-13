function solutions = all_solution(algo,pro,m,d,runNo)
    numResult = 1;    % 取り出したい result のインデックス
    solutions.decs  = [];  % 決定変数を格納する配列
    solutions.objs = [];  % 目的関数を格納する配列
    solutions.cons = [];  % 拘束違反量を格納する配列

    solutions.initDecs  = []; 
    solutions.initObjs = [];
    solutions.initCons = [];

    initPopSize = 100;
    for a = 1:length(algo)
        algorithm = algo{a};
        for i = 1:runNo
            data   = load(sprintf('Data/NSGAII_all/%s/%s_%s_M%d_D%d_%d.mat',algorithm,algorithm,pro,m,d, i),'result');
            result = data.result;
            
            popData = result{numResult,2};  
            solutions.decs  = [solutions.decs;  popData.decs];
            solutions.objs = [solutions.objs; popData.objs];
            solutions.cons = [solutions.cons; popData.cons];
            dec = [popData.decs];
            obj = [popData.objs];
            con = [popData.cons];
            selIdx = 1:initPopSize;
            solutions.initDecs  = [solutions.initDecs;  dec(selIdx,:)];
            solutions.initObjs = [solutions.initObjs; obj(selIdx,:)];
            solutions.initCons = [solutions.initCons; con(selIdx,:)];
        end
    end
    disp(solutions)
    
    if ~isempty(solutions.decs)
        [uniqueDecs, idxUnique] = unique(solutions.decs, 'rows');
        solutions.decs  = uniqueDecs;
        solutions.objs = solutions.objs(idxUnique, :);
        solutions.cons = solutions.cons(idxUnique, :);
    end

    % --- 重複解を除外（初期個体群） ---
    if ~isempty(solutions.initDecs)
        [uniqueInitDecs, idxInitUnique] = unique(solutions.initDecs, 'rows');
        solutions.initDecs  = uniqueInitDecs;
        solutions.initObjs = solutions.initObjs(idxInitUnique, :);
        solutions.initCons = solutions.initCons(idxInitUnique, :);
    end
    disp(solutions)
    %{
    feasible_index   = find(sum(max(0, solutions.cons), 2) == 0);
    infeasible_index = setdiff(1:size(solutions.objs, 1), feasible_index);

    figure;
    scatter(solutions.objs(infeasible_index, 1), solutions.objs(infeasible_index, 2), ...
        10, [0.8 0.8 1], 'filled', 'MarkerEdgeColor', 'blue');
    hold on;
    scatter(solutions.objs(feasible_index, 1), solutions.objs(feasible_index, 2), ...
        10, [1 0.8 0.8], 'filled', 'MarkerEdgeColor', 'red');
    hold off;
    %}
end