function [Population,SubPopulation,SubPop_Constraint_value,SubPop_Scalarizing_value,Del,W] = updateWeight_dynamic(Population,SubPopulation,SubPop_Constraint_value,SubPop_Scalarizing_value,W,Z,EP,change_ind)
% Delete overcrowded subproblems and add new subproblems

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [~,M] = size(Population.objs);
    [N,S] = size(SubPopulation);
    H = 0.25;
    L = 0.05;
    nus = length(change_ind)/N;
    nus = (2*(H-L)*nus-H+2*L)*N;
 %% Update the current population by EP
    feasible_index = find(sum(max(0,Population.cons),2) == 0);
    Pop_fea = Population(feasible_index);
    Combine = [Pop_fea,EP];
    CombineObj = abs(Combine.objs-repmat(Z,length(Combine),1));
    g = zeros(length(Combine),size(W,1));
    for i = 1 : size(W,1)
        g(:,i) = max(CombineObj.*repmat(W(i,:),length(Combine),1),[],2);
    end
    % Choose the best solution for each subproblem
    [~,best]   = min(g,[],1);
    Population = Combine(best);
    % Choose the best solution for each subproblem
 
    %% Delete the overcrowded subproblems
    Dis = pdist2(Population.objs,Population.objs);
    Dis(logical(eye(length(Dis)))) = inf;
    Del = false(1,length(Population));
    t = 0;
    while t < min(nus,length(EP))
        Remain = find(~Del);
        subDis = sort(Dis(Remain,Remain),2);
        [~,worst] = sort(prod(subDis(:,1:min(M,length(Remain))),2));
        i = 1;
        while Remain(worst(i)) ~= change_ind
            i = i+1;
            t = t+1;
        end
        Del(Remain(worst(i))) = true;
        t = t+1;
    end
    Population = Population(~Del);
    W = W(~Del,:);
    SubPopulation = SubPopulation(~Del,:);
    SubPop_Constraint_value = SubPop_Constraint_value(~Del,:);
    SubPop_Scalarizing_value = SubPop_Scalarizing_value(~Del,:);
    
    %% Add new subproblems
    % Determine the new solutions be added
    Combine  = [Population,EP];
    Selected = false(1,length(Combine));
    Selected(1:length(Population)) = true;
    Dis = pdist2(Combine.objs,Combine.objs);
    Dis(logical(eye(length(Dis)))) = inf;
    
    while sum(Selected) < min(N,length(Selected))
        subDis = sort(Dis(~Selected,Selected),2);
        [~,best] = max(prod(subDis(:,1:min(M,size(subDis,2))),2));
        Remain = find(~Selected);
        Selected(Remain(best)) = true;
    end
    
    % Add new subproblems
    newObjs = EP(Selected(length(Population)+1:end)).objs;
    temp    = 1./(newObjs-repmat(Z,size(newObjs,1),1));
    W = [W;temp./repmat(sum(temp,2),1,size(temp,2))];
    % Add new solutions
    Population = Combine(Selected);
end