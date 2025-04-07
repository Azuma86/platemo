%=======================================================================
% 「data」フォルダ内のアルゴリズム名のフォルダ (例: 'NSGAIISR') のパス
algFolderPath = fullfile('data','NSGAIICVM');  

% 特定の文字列を含むファイルのみを対象にする (例: 'RWMOP11')
searchString = 'RWMOP50';

% algFolderPath に含まれる .mat ファイルリストを取得
matList = dir(fullfile(algFolderPath, '*.mat'));

% 取得した全 MAT ファイルに対して処理を行う
for i = 1:numel(matList)
    % ファイル名を取得
    matFileName = matList(i).name;

    % ファイル名に指定文字列 searchString が含まれるかチェック
    if contains(matFileName, searchString)

        % フルパスを作成
        matFilePath = fullfile(algFolderPath, matFileName);

        % MATファイルの中身を読み込む
        fileData = load(matFilePath);

        % metric という変数が含まれている場合，削除して上書き保存
        if isfield(fileData, 'metric')
            fileData = rmfield(fileData, 'metric');   % 構造体からフィールドを削除
            save(matFilePath, '-struct', 'fileData'); % 上書き保存
            fprintf('Deleted "metric" from %s\n', matFileName);
        end
    end
end