classdef PM_one < ALGORITHM
    % <multi> <real/integer/label/binary/permutation> <constrained/none>
    methods
        function main(Algorithm,Problem)
            pn = Problem.number;
            rng(125)
            N = Problem.N;
            D = Problem.D;

            % 制約インデックスの指定
            targetConstraints = [4]; 
            
            % LHS初期化
            lhsSample = lhsdesign(N, D);
            lowerBounds = repmat(Problem.lower, N, 1);
            upperBounds = repmat(Problem.upper, N, 1);
            PopDec = lowerBounds + lhsSample .* (upperBounds - lowerBounds);
            Population = Problem.Evaluation(PopDec);

            % 突然変異関連パラメータ
            proM = 0.9;
            disM = 20;
            
            PrevPopulation = Population;
            stableCount = zeros(1,N);
            stableThreshold = 20;

            % ===============================
            % CSVファイルの初期化
            % ===============================
            csvFilename = sprintf('/Users/azumayuki/Documents/LONs/data09-20-pre/RWMOP%d_data.csv',pn);
            objCount = length(Population(1).objs);
            conCount = length(Population(1).cons);
            decHeader = arrayfun(@(x) sprintf('X_%d',x), 1:D, 'UniformOutput', false);
            objHeader = arrayfun(@(x) sprintf('Obj_%d',x), 1:objCount, 'UniformOutput', false);
            conHeader = arrayfun(@(x) sprintf('Con_%d',x), 1:conCount, 'UniformOutput', false);
            header = [{'Gen', 'ID'}, decHeader{:}, objHeader{:}, conHeader{:}];
            
            fid = fopen(csvFilename, 'w');
            fprintf(fid, '%s,', header{1:end-1});
            fprintf(fid, '%s\n', header{end});
            fclose(fid);

            currentGeneration = 0;

            % ===============================
            % メインループ
            % ===============================
            while Algorithm.NotTerminated(Population)
                currentGeneration = currentGeneration + 1;

                % 現世代個体をCSVへ追記
                fid = fopen(csvFilename, 'a');
                for i = 1:length(Population)
                    individual = Population(i);
                    dec = individual.decs;
                    x = [];
                    for t = 1:length(x)
                        dec(x(t)) = round(dec(x(t)));
                    end
                    decStr = strjoin(arrayfun(@num2str, dec, 'UniformOutput', false), ',');
                    objStr = strjoin(arrayfun(@num2str, individual.objs, 'UniformOutput', false), ',');
                    conStr = strjoin(arrayfun(@num2str, individual.cons, 'UniformOutput', false), ',');
                    fprintf(fid, '%d,%d,%s,%s,%s\n', currentGeneration, i, decStr, objStr, conStr);
                end
                fclose(fid);

                % ===============================
                % 次世代生成：1+10個体から最良を選択
                % ===============================
                Offspring = Population;
                for i = 1:N
                    Parent = Population(i);
                    if stableCount(i) < stableThreshold
                        % 親個体 + 10子 → 計11個体
                        Candidates = repmat(Parent,1,11);
                        for j = 1:10
                            mutatedDec = polynomialMutation(Parent.decs, Problem.lower, Problem.upper, proM, disM);
                            Candidates(j+1) = Problem.Evaluation(mutatedDec);
                        end
                        
                        % ========= 修正ポイント: 指定制約だけを足し合わせる =========
                        % Candidates(k).cons は (1×conCount)
                        % -> cat(1, ...) で (11×conCount) の行列を生成
                        PopCon = cat(1, Candidates.cons);  % (11, conCount)
                        
                        % targetConstraints=2 なら PopCon(:,2) が制約#2
                        selectedCons = PopCon(:, targetConstraints);  % (11, 1)
                        
                        % 正の部分だけ足し合わせ → CV_values(11×1)
                        CV_values = sum(max(0, selectedCons),2);
                        
                        % もし複数制約を加算したい場合、sum(...,2)などで行方向に和
                        % ここでは1本だけなのでそのまま
                        % 例) CV_values = sum(max(0, selectedCons), 2);
                        
                        [~, bestIdx] = min(CV_values);
                        Offspring(i) = Candidates(bestIdx);
                    else
                        Offspring(i) = Parent;
                    end
                end

                Population = Offspring;

                % 安定度カウンタ更新
                for i = 1:N
                    if isequal(Population(i).decs, PrevPopulation(i).decs)
                        stableCount(i) = stableCount(i) + 1;
                    else
                        stableCount(i) = 0;
                    end
                end
                if all(stableCount >= stableThreshold)
                    break;
                end
                PrevPopulation = Population;
            end
        end
    end
end

function mutatedOffspring = polynomialMutation(Offspring, lower, upper, proM, disM)
    D = length(Offspring);
    for i = 1:D
        if rand < proM
            mu = rand;
            yl = lower(i);
            yu = upper(i);
            if mu <= 0.5
                deltaq = ((2*mu) + (1 - 2*mu)*((1-(Offspring(i)-yl)/(yu-yl))^(disM+1)))^(1/(disM+1)) - 1;
                Offspring(i) = Offspring(i) + deltaq*(yu-yl);
            else
                deltaq = 1 - ((2*(1-mu)) + (2*(mu-0.5)*((1-(yu-Offspring(i))/(yu-yl))^(disM+1))))^(1/(disM+1));
                Offspring(i) = Offspring(i) + deltaq*(yu-yl);
            end
        end
    end
    mutatedOffspring = min(max(Offspring, lower), upper);
end