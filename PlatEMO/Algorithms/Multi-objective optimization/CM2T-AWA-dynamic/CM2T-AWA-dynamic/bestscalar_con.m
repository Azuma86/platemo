function r = bestscalar_con(Scalarzing_value,Constraint_value,current_best,gen,last_gen)
    [Minimize_Constraint_value,~] = min(Constraint_value, [], 2);
    bestscalar = zeros(1,length(Scalarzing_value));
    for i = 1 : length(Scalarzing_value)
        minimize_constraint_index_insubpop = find(Constraint_value(i,:) == Minimize_Constraint_value(i));
        if length(minimize_constraint_index_insubpop) > 1
            [~,index] = min(Scalarzing_value(i,minimize_constraint_index_insubpop));
            bestscalar(i) = Scalarzing_value(i,minimize_constraint_index_insubpop(index(1)));
        else
            bestscalar(i) = Scalarzing_value(i,minimize_constraint_index_insubpop);
        end
    end
    % 現在の世代からlast_gen世代前までのスカラー化関数値を保存
    if gen <= last_gen
        current_best(gen, :) = bestscalar;
    else
        % 一番古い世代を置き換える
        current_best = circshift(current_best, -1, 1); % 上方向にシフト
        current_best(last_gen, :) = bestscalar; % 新しい世代を追加
    end
    r = current_best;
end