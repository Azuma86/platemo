function r = bestscalar(Scalarzing_value,current_best,gen,last_gen)
    [bestscalar,~] = min(Scalarzing_value, [], 2);
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