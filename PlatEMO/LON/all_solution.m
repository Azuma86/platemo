function solutions = all_solution(algo,pro,m,d,runNo)
    numResult = 1;    % 取り出したい result のインデックス
    solutions.dec  = [];  % 決定変数を格納する配列
    solutions.objs = [];  % 目的関数を格納する配列
    solutions.cons = [];  % 拘束違反量を格納する配列
    for a = 1:length(algo)
        algorithm = algo{a};
        for i = 1:runNo
            data   = load(sprintf('%s_%s_M%d_D%d_%d.mat',algorithm,pro,m,d, i),'result');
            result = data.result;
            
            popData = result{numResult,2};  
            solutions.dec  = [solutions.dec;  popData.decs];
            solutions.objs = [solutions.objs; popData.objs];
            solutions.cons = [solutions.cons; popData.cons];
        end
    end
    disp(solutions)
    
    if ~isempty(solutions.dec)
        [uniqueDec, idxUnique] = unique(solutions.dec, 'rows');
        solutions.dec  = uniqueDec;
        solutions.objs = solutions.objs(idxUnique, :);
        solutions.cons = solutions.cons(idxUnique, :);
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