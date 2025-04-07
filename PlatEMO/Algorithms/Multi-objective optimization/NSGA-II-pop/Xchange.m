function r = Xchange(x,current_X,gen,last_gen)
    if gen <= last_gen
        current_X(gen, :) = x;
    else
        % 一番古い世代を置き換える
        current_X = circshift(current_X, -1, 1); % 上方向にシフト
        current_X(last_gen, :) = x; % 新しい世代を追加
    end
    r = current_X;
end