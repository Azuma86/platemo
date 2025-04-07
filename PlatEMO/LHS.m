
    % problemIndex: RWMOP問題番号 (1～50など)
    % N: サンプリング数
    % filenameX: 決定変数Xを書き出すCSVファイル名
    % filenameC: 制約違反量Consを書き出すCSVファイル名

    % --- 1) 問題をインスタンス化 ---
    % 例）"RWMOP15" を呼び出す場合は下記のようにevalで文字列生成
    problemIndex = 3;
    N = 100000;
    filenameX = sprintf('/Users/azumayuki/Documents/LONs/data_LHS/X_RWMOP%d.csv',problemIndex);
    filenameC = sprintf('/Users/azumayuki/Documents/LONs/data_LHS/Cons_RWMOP%d.csv',problemIndex);
    eval(['Problem = RWMOP', num2str(problemIndex), '();']);
    Problem.Setting();  % ここでProblem.lower, Problem.upperなどが決定

    D = Problem.D;       % 次元数
    lb = Problem.lower;  % 下限ベクトル
    ub = Problem.upper;  % 上限ベクトル

    % --- 2) LHSサンプリング ---
    % lhsdesign(N, D): [0,1]の範囲でLHSを行った後，下限～上限にスケーリング
    latinCube = lhsdesign(N, D);
    X = bsxfun(@plus, lb, bsxfun(@times, latinCube, (ub - lb)));

    % --- 3) 評価 ---
    Population = Problem.Evaluation(X);
    Cons = Population.cons;   % 制約違反量 (N×C)
    CV = sum(max(0,Cons),1);
    disp(length(find(CV == 0)))
    % --- 4) ファイル出力 ---
    % 行列をCSVに書き出し
    % （行列サイズが大きい場合はwritematrixやdlmwriteの方が高速なことも）
    %writematrix(X, filenameX);
    %writematrix(Cons, filenameC);

    %fprintf('Done sampling RWMOP%d: saved %d samples.\n', problemIndex, N);