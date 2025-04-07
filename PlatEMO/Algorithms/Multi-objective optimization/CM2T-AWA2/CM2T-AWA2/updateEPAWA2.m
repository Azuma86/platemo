function EP = updateEPAWA2(EP,Offsprings,nEP)
% Update the external population

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Select the non-dominated solutions
    EP = [EP,Offsprings];
    feasible_index = find(sum(max(0,EP.cons),2) == 0);
    EP = EP(feasible_index);
    EP = EP(NDSort(EP.objs,1)==1);
    [N,M] = size(EP.objs);
    %% Delete the overcrowded solutions
    Dis = pdist2(EP.objs,EP.objs);
    Dis(logical(eye(length(Dis)))) = inf;
    Del = false(1,N);
    while sum(Del) < N-nEP
        Remain = find(~Del);
        density = estimateDensity(Remain,Dis(Remain,Remain),M);
        [~,worst] = max(density);
        Del(Remain(worst)) = true;
    end
    EP = EP(~Del);
end

function density = estimateDensity(Population,Dis,M)
    % 各解の密度を計算
    density = zeros(1, length(Population));
    for i = 1:length(Population)
        sortedDistances = sort(Dis(i, :));
        density(i) = sum(1 ./ (sortedDistances(1:min(M,length(Population)))+1e-6)); % 自分自身を除外
    end
end