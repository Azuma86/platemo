classdef NSGAIISR < ALGORITHM
    % <multi> <real/integer/label/binary/permutation> <constrained/none>
    % NSGA-II with Infeasible Elitists and Stochastic Ranking Selection (IE-SRS)
    %
    % Modified Step1:
    %  - Non-dominated solutions are identified among ALL solutions (feasible + infeasible)
    %  - Then only the feasible ones (CV=0) among them are kept.
    %  - If these exceed N, use crowding distance truncation and skip Step3~.
    properties
        Pf = 0.45;  % Probability threshold for objective-based comparison (SR)
        R  = SOLUTION.empty(); % Infeasible elitists set
    end
    
    methods
        function main(Algorithm, Problem)
            pn = Problem.number;
            %% Initialization
            Population = Problem.Initialization();        
            Algorithm.R = SOLUTION.empty(); 
            f = 0;
            gen  = 1;
            [Population, FrontNo, CrowdDis] = EnvironmentalSelectionIE_SRS(...
                Population, Algorithm.R, Problem.N, Algorithm.Pf);
            filename = sprintf('/Users/azumayuki/Documents/LONs/feasible_ratio_data/RWMOP%d_SR.csv', pn);
            filename2 = sprintf('/Users/azumayuki/Documents/LONs/RWMOP%d_SR_first.csv', pn);

            f = logFeasibleRatio(Population, gen, filename, f, filename2);
            arch = updateEP2(Population);
            %% Optimization loop
            while Algorithm.NotTerminated(Population)
                gen = gen + 1;
                % (1) Generate mating pool by binary tournament
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, -CrowdDis);
                % (2) Variation operators
                Offspring  = OperatorGA(Problem, Population(MatingPool));
                arch = updateEP2([arch,Offspring]);
                % (3) Update infeasible elitists set
                Algorithm.R = UpdateInfeasibleElitists(Algorithm.R, Population, Problem.N);
                
                % (4) Environmental Selection with IE-SRS
                [Population, FrontNo, CrowdDis] = EnvironmentalSelectionIE_SRS(...
                    [Population, Offspring, Algorithm.R], Algorithm.R, Problem.N, Algorithm.Pf);
                f = logFeasibleRatio(Population, gen, filename, f, filename2);
                if Problem.FE >= Problem.maxFE
                   arch = archive(arch,Problem.N*2);
                   Population = arch;
                end
            end
        end
    end
end

%% ========================================================================= %%
function [PopNext, FrontNo, CrowdDis] = EnvironmentalSelectionIE_SRS(...
                                        PopAll, R, N, Pf)
% ENVIRONMENTALSELECTIONIE_SRS
%
% Step1: Among ALL solutions (PopAll), identify non-dominated solutions => pick only feasible ones
% Step2: If they exceed N, crowding distance truncation => done (skip step3~)
% Step3: Otherwise, define the rest (dominated feasible + all infeasible), do SR-based sorting, etc.
%

%% Step1: Identify "non-dominated solutions among ALL" and pick only feasible
CV_all   = sum(max(0, PopAll.cons), 2);
[FrontNoAll, ~]  = NDSort(PopAll.objs, length(PopAll));   % Non-dominated sort (all)
NonDomIndex = (FrontNoAll == 1);                          % Indices of front=1
PopNonDom   = PopAll(NonDomIndex);                        % All non-dominated (feasible+infeasible)

% そこから可行解(CV=0)のみを抽出
CV_nonDom       = CV_all(NonDomIndex);
feasibleInND    = (CV_nonDom == 0);
PopFeasibleND   = PopNonDom(feasibleInND);

% これらがStep1で選ばれた個体
PopNext = PopFeasibleND;

%% Step2: If Step1 solutions exceed N => crowding distance truncation & finish
if length(PopNext) > N
    [Ftemp, ~]  = NDSort(PopNext.objs, length(PopNext));
    Ctemp       = CrowdingDistance(PopNext.objs, Ftemp);
    exceed      = length(PopNext) - N;
    % Remove the smallest crowdDist first
    [~, so]     = sort(Ctemp, 'ascend');
    PopNext(so(1:exceed)) = [];
    
    % ここで次世代がN個確定 => 以降のStep3は行う意味がない
    [FrontNo, ~] = NDSort(PopNext.objs, N);
    CrowdDis     = CrowdingDistance(PopNext.objs, FrontNo);
    return;  % <<<========  Step3以降はスキップして終了
end

%% Step3: ここに来た段階で、まだPopNext.size <= N
%         => 残り枠を埋めるために、dominated feasible + infeasible をSRで評価
IndexInNext = false(1,length(PopAll));
for i = 1 : length(PopNext)
    idx = find(PopAll == PopNext(i),1);
    if ~isempty(idx)
        IndexInNext(idx) = true;
    end
end
Pr = PopAll(~IndexInNext);  % the "rest"

% --------- 以下、SR-basedバブルソートで Pr を並び替える & 取り込み --------- %

%% (a) Prepare feasible-solution ranks/crowding inside Pr
Pr_CV = sum(max(0, Pr.cons), 2);

Feas_idx   = find(Pr_CV == 0);
Objs_feas  = cat(1, Pr(Feas_idx).objs);

if ~isempty(Feas_idx)
    [Fno_loc, ~] = NDSort(Objs_feas, length(Feas_idx));
    Cdis_loc     = CrowdingDistance(Objs_feas, Fno_loc);
else
    Fno_loc  = [];
    Cdis_loc = [];
end

% Assign rank/crowd to each individual in Pr
rankMap  = inf(1,length(Pr));
crowdMap = zeros(1,length(Pr));
for k = 1:length(Feas_idx)
    idxSol = Feas_idx(k);
    rankMap(idxSol)  = Fno_loc(k);
    crowdMap(idxSol) = Cdis_loc(k);
end

% (b) SR-based bubble sort with multiple passes
maxPass = 2 * length(Pr);
K = 1;
sflag = 1;
while (K <= maxPass) && (sflag == 1)
    sflag = 0;
    i = 1;
    while i < length(Pr)
        Pi   = Pr(i);
        Pi1  = Pr(i+1);
        CVi  = Pr_CV(i);
        CVi1 = Pr_CV(i+1);
        
        if CVi == 0 && CVi1 == 0
            % Both feasible => compare by (rank, -crowding)
            if ~CompareByRankCrowd(rankMap(i), crowdMap(i), ...
                                   rankMap(i+1), crowdMap(i+1))
                [Pr(i), Pr(i+1)] = deal(Pi1, Pi);
                [Pr_CV(i), Pr_CV(i+1)] = deal(CVi1, CVi);
                [rankMap(i), rankMap(i+1)]   = deal(rankMap(i+1), rankMap(i));
                [crowdMap(i), crowdMap(i+1)] = deal(crowdMap(i+1), crowdMap(i));
                sflag = 1;
            end
        else
            % At least one is infeasible => Stochastic Ranking
            mu = rand();
            if mu <= Pf
                % objective-based => dominance
                if ~dominates(Pi.objs, Pi1.objs)
                    [Pr(i), Pr(i+1)] = deal(Pi1, Pi);
                    [Pr_CV(i), Pr_CV(i+1)] = deal(CVi1, CVi);
                    [rankMap(i), rankMap(i+1)]   = deal(rankMap(i+1), rankMap(i));
                    [crowdMap(i), crowdMap(i+1)] = deal(crowdMap(i+1), crowdMap(i));
                    sflag = 1;
                end
            else
                % constraint-based => CV
                if CVi > CVi1
                    [Pr(i), Pr(i+1)] = deal(Pi1, Pi);
                    [Pr_CV(i), Pr_CV(i+1)] = deal(CVi1, CVi);
                    [rankMap(i), rankMap(i+1)]   = deal(rankMap(i+1), rankMap(i));
                    [crowdMap(i), crowdMap(i+1)] = deal(crowdMap(i+1), crowdMap(i));
                    sflag = 1;
                end
            end
        end
        
        i = i + 1;
    end
    K = K + 1;
end

%% (c) Insert top from Pr until we fill up PopNext
remainSlots = N - length(PopNext);
if remainSlots > 0
    toAdd = Pr(1 : min(remainSlots, length(Pr)));
    PopNext = [PopNext, toAdd];
end

%% (d) If exceed => crowding distance truncation
if length(PopNext) > N
    [Ftmp, ~]  = NDSort(PopNext.objs, length(PopNext));
    Ctmp       = CrowdingDistance(PopNext.objs, Ftmp);
    exceed = length(PopNext) - N;
    [~, so] = sort(Ctmp, 'ascend');
    PopNext(so(1:exceed)) = [];
end

%% Final: compute FrontNo, CrowdDis for the returned population
[FrontNo, ~] = NDSort(PopNext.objs, N);
CrowdDis     = CrowdingDistance(PopNext.objs, FrontNo);

end

%% ========================================================================= %%
function R = UpdateInfeasibleElitists(R, Population, MaxSize)
% Keep infeasible solutions in R if they are not dominated by existing R
for i = 1:length(Population)
    cv = sum(max(0, Population(i).cons));
    if cv > 0
        dominated = false;
        for j = 1:length(R)
            if dominates(R(j).objs, Population(i).objs) && ...
                    (sum(max(0, R(j).cons)) <= cv)
                dominated = true;
                break;
            end
        end
        if ~dominated
            R = [R, Population(i)];
            if length(R) > MaxSize
                % Remove the worst CV
                CV_all = sum(max(0, cat(1,R.cons)),2);
                [~, idxWorst] = max(CV_all);
                R(idxWorst) = [];
            end
        end
    end
end
end

%% ========================================================================= %%
function better = CompareByRankCrowd(rank1, crowd1, rank2, crowd2)
% Compare feasible solutions by (rank, crowding)
if rank1 < rank2
    better = true;
elseif rank1 > rank2
    better = false;
else
    % same rank => prefer larger crowding
    better = (crowd1 > crowd2);
end
end

%% ========================================================================= %%
function flag = dominates(obj1, obj2)
% Pareto dominance check
flag = all(obj1 <= obj2) && any(obj1 < obj2);
end

function t = logFeasibleRatio(Pop, gen, filename,f ,filename2)
    CV = arrayfun(@(p) sum(max(0,p.cons)), Pop);
    feasibleCount = sum(CV == 0);
    
    ratio = feasibleCount / length(Pop);
    if f == 0
        if ratio > 0
            f = 1;
            fid = fopen(filename2,'a');
            fprintf(fid,"%d,%f\n",gen);
            fclose(fid);
        end
    end
    t = f;
    fid = fopen(filename,'a');
    fprintf(fid,"%d,%f\n",gen,ratio);
    fclose(fid);
end