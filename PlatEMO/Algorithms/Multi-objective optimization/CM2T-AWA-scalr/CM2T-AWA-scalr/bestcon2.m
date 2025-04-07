function r = bestcon2(Constraint_value,current_best,gen,last_gen)
    [Minimize_Constraint_value,~] = min(Constraint_value, [], 2);
    % 現在の世代からlast_gen世代前までのスカラー化関数値を保存
    if gen <= last_gen
        current_best(gen, :) = Minimize_Constraint_value;
    end
    r = current_best;
end