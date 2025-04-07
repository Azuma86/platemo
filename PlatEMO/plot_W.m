% Load weight vectors and reference point
W = readmatrix('w_vectors.txt');  % 重みベクトル
Z = readmatrix('z_point.txt');    % 参照点（理想点）

% 2D plot
if size(W, 2) == 2
    hold on;
    extend_length = 5; % 突き抜ける長さのスケール

    % 各ベクトルを突き抜ける形で描画
    for i = 1:size(W, 1)
        % 重みベクトルの方向を単位ベクトル化
        direction = W(i,:) - Z;
        direction = direction / norm(direction);

        % 直線の始点と終点
        start_point = Z;
        end_point = W(i,:) + extend_length * direction; % 突き抜ける

        % 線を描画
        plot([start_point(1), end_point(1)], [start_point(2), end_point(2)], 'r-', 'LineWidth' ,0.5);
    end
    hold off;
end