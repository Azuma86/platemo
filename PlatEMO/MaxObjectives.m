function [minObj, maxObj] = MaxObjectives()

    %% ==== パラメータ設定 ====
    %algList = {"NSGAIICDP","NSGAIISP","NSGAIISR","NSGAIIMo","NSGAIIeps"};  % 解析したいアルゴリズム名リスト
    algList = {"NSGAIICDP","NSGAIIcount","NSGAIIchev","NSGAIICVM"};  % 解析したいアルゴリズム名リスト
    Pro     = 50;  % 問題番号 (RWMOPなど)
    runNo   = 31;  % 試行回数
    M       = 2;   % 目的関数の次元
    D       = 6;  % 変数次元 (必要なら変更)
    numResult = 1; % result{row, :} で取り出す行番号 (PlatEMOのresultの保存方法に合わせて設定)

    %% ==== すべての解を格納する配列 ====
    allObj = [];  % ここに全アルゴリズム・全試行の全解の目的値を集める

    %% ==== 全アルゴリズム・全試行をループして解を取得 ====
    for alg = 1:length(algList)
        algName = algList{alg};

        for i = 1:runNo
            fileName = sprintf('%s_RWMOP%d_M%d_D%d_%d.mat', algName, Pro, M, D, i);

            % ファイルが存在する場合のみ読み込み
            if exist(fileName, 'file')
                data = load(fileName, 'result');  % 'result'のみ読み込み

                if isfield(data,'result')
                    thisResult = data.result;
                    % numResult行が存在しているかチェック
                    if size(thisResult,1) >= numResult
                        % 個体群 (Population など) を取り出す想定
                        Population = thisResult{numResult, 2};
                        Obj = Population.objs; % 目的関数値行列(N×M)

                        % NaN でないかチェック (NaNが混在していると max/min でエラーになる可能性)
                        if ~any(isnan(Obj(:)))
                            allObj = [allObj; Obj]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end

    %% ==== 非劣解の抽出 ====
    % ここでは PlatEMO の NDSort を利用できる想定で書いています
    if ~isempty(allObj)
        ND = NDSort(allObj, 1) == 1;       % 非劣解(ND解)のインデックス (true/false)
        ndObj = allObj(ND, :);            % 非劣解の目的関数値

        %% ==== 最小値 & 最大値を求める ====
        minObj = min(ndObj,[],1);
        maxObj = max(ndObj,[],1);
    else
        % 1つも解が読み込めなかった場合の対処: ここでは空で返す
        warning('解が見つかりませんでした．');
        minObj = [];
        maxObj = [];
    end

    %% ==== 結果の表示 (必要に応じて) ====
    disp('===== 非劣解から求めた各目的関数の最小値・最大値 =====');
    for m = 1:length(minObj)
        fprintf('Objective %d: min = %.4f, max = %.4f\n', ...
            m, minObj(m), maxObj(m));
    end
    fprintf('R = [%8f %8f;%8f %8f];\n',minObj(1),minObj(2),maxObj(1),maxObj(2))
    %fprintf('[%8f %8f %8f;%8f %8f %8f]\n',minObj(1),minObj(2),minObj(3),maxObj(1),maxObj(2),maxObj(3))
    %fprintf('[%8f %8f %8f %8f %8f;%8f %8f %8f %8f %8f]\n',minObj(1),minObj(2),minObj(3),minObj(4),minObj(5),maxObj(1),maxObj(2),maxObj(3),maxObj(4),maxObj(5))
end