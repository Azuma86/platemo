function neighbor_search(solutions,pro)
   
    k        = 3;      % k近傍
    numStart = 1000;     % ランダムに選ぶ開始点数    

    N = size(solutions.dec,1);
    D = size(solutions.dec,2);
    C = size(solutions.cons,2);
    
    % 制約違反総量
    cvAll = sum(max(0, solutions.cons), 2);  % N×1
    
    % 近傍探索用の k-d-tree モデルを作成
    % NSMethod='kdtree' はデフォルトで2～10次元くらいの時に有効
    % Dが大きいと性能に注意
    Mdl = createns(solutions.dec,'NSMethod','kdtree');
    
    % 出発点をランダムに選択
    rng(0,'twister');
    startIDs = randperm(N, min(numStart,N));
    
    % ========== 結果格納用 (テーブルに出す) ===========
    % Python 側での可視化に合わせた列名を用意
    varNames = ["ID","Gen", ...
                arrayfun(@(d) sprintf('X_%d', d), 1:D, 'UniformOutput',false), ...
                arrayfun(@(c) sprintf('Con_%d',c), 1:C, 'UniformOutput',false)];
    %varNames = [varNames{:}];  % cellstr 化
    
    allRecords = [];  % 後で table に変換
    
    % ========== 局所探索 (各開始点で繰り返し) ===========
    for s = 1:length(startIDs)
        currentID  = startIDs(s);   % 現在の点のインデックス
        cvCurrent  = cvAll(currentID);
        decCurrent = solutions.dec(currentID,:);
        consCurrent= solutions.cons(currentID,:);
        
        gen = 0;
        tempRecords = [];
        
        while true
            % 1) 記録 (ID, Gen, X_..., Con_...)
            rowData = [s, gen, decCurrent, consCurrent];
            tempRecords = [tempRecords; rowData];
            %
            % 2) k近傍探索 (knnsearch)
            %    k+1 を取得: 1つ目が自分自身, 残りk個が近傍
            [neighborIDs, distVals] = knnsearch(Mdl, decCurrent, 'K', k+1);
            
            neighborIDs(1) = [];
            distVals(1) = [];
            % 3) 近傍の CV をチェックし，現在より小さなものを探す
            cvNeighbors = cvAll(neighborIDs);
            %}

            %{
            distVec = pdist2(decCur, solutions.dec);  % 1×N
            [distSort, idxSort] = sort(distVec, 'ascend');
            
            % idxSort(1) は自分自身(距離0)
            neighborIDs = idxSort(2 : k+1);  % 最も近い k 個のインデックス
            
            % 3) 近傍の中で「CVが現在より小さい点」を探す
            cvNeighbors = cvAll(neighborIDs);
            %}

            betterIdx   = neighborIDs(cvNeighbors < cvCurrent);
            
            if isempty(betterIdx)
                % 近傍により良い(CVが低い)解がなければ停止
                break;
            end
            
            % 4) その中で最もCVが小さいものを選ぶ
            [~, idxMin] = min(cvAll(betterIdx));
            bestID = betterIdx(idxMin);
            
            % 移動
            currentID  = bestID;
            cvCurrent  = cvAll(bestID);
            decCurrent = solutions.dec(bestID,:);
            consCurrent= solutions.cons(bestID,:);
            
            gen = gen + 1;
        end
        
        % この開始点の探索経路をまとめて allRecords に追加
        allRecords = [allRecords; tempRecords];
    end
    % ========== テーブル化 & CSV 出力 ===========
    T = array2table(allRecords, 'VariableNames', varNames);
    writetable(T, sprintf('/Users/azumayuki/Documents/LON/data09-20-pre/local_search%s.csv',pro));
    fprintf('Local search finished. %d records saved.\n', size(T,1));
end