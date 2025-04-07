% LIRCMOP7クラスを使ってパレートフロントをプロットするスクリプト

% LIRCMOP7クラスのインスタンスを生成
problem = MW7();
problem.Setting(); % 問題の設定

% パレートフロントの点を取得
N = 1000; % プロットする点の数
R = problem.GetOptimum(N); % 最適解を取得
if length(R(1,:)) == 2
    % パレートフロントをプロット
    scatter(R(:,1), R(:,2), 'filled', 'MarkerFaceColor', 'magenta');
else
    scatter3(R(:,1),R(:,2),R(:,3), 'filled', 'MarkerFaceColor', 'red');
end
