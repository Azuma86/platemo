function r = ideal(ideal_point,current_ideal,gen,last_gen)
    % 現在の世代からlast_gen世代前までのスカラー化関数値を保存
    if gen <= last_gen
        current_ideal(gen, :) = ideal_point;
    else
        % 一番古い世代を置き換える
        current_ideal = circshift(current_ideal, -1, 1); % 上方向にシフト
        current_ideal(last_gen, :) = ideal_point; % 新しい世代を追加
    end
    r = current_ideal;
end