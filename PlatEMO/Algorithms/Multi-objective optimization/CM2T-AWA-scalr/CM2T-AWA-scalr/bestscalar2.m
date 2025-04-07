function r = bestscalar2(Scalarzing_value,current_best,gen,last_gen)
    [bestscalar,~] = min(Scalarzing_value, [], 2);
    % 現在の世代からlast_gen世代前までのスカラー化関数値を保存
    if gen <= last_gen
        current_best(gen, :) = bestscalar;
    end
    r = current_best;
end