function [Del,W] = updateWeightAWA(Population,SubPopulation,W,Z,EP,nus,T)
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

    %% Update the current population by EP
    % Calculate the function value of each solution in Population or EP on
    % each subproblem in W
    Combine = [Population,EP];
    CombineObj = abs(Combine.objs-repmat(Z,length(Combine),1));
    g = zeros(length(Combine),size(W,1));
    for i = 1 : size(W,1)
        g(:,i) = max(CombineObj.*repmat(W(i,:),length(Combine),1),[],2);
    end
    % Choose the best solution for each subproblem
    [~,best]   = min(g,[],1);
    Population = Combine(best);
    
    SubPopulation_vec = reshape(SubPopulation',1,[]);
    
    %% Delete the overcrowded subproblems
    Dis = pdist2(SubPopulation_vec.objs,SubPopulation_vec.objs);
    %Dis(logical(eye(length(Dis)))) = inf;
    Dis = sort(Dis,2);
    Dis = Dis(:,1:S*T);
    Dis = sum(Dis,2);
    Dis = reshape(Dis,S,N);
    Dis = sum(Dis);
    Del = false(1,length(Population));
    while sum(Del) < min(nus,length(EP))
        Remain = find(~Del);
        [~,minIndex] = min(Dis(Remain));
        %subDis = sort(Dis(Remain,Remain),2);
        %min_value = min(prod(subDis(:,1:min(M,length(Remain))),2));
        %min_indices = find(prod(subDis(:,1:min(M,length(Remain))),2) == min_value);
        %random_index = min_indices(randi(length(min_indices)));
        Del(Remain(minIndex)) = true;
    end
    Population = Population(~Del);
    W = W(~Del,:);
    
    %% Add new subproblems
    % Determine the new solutions be added
    Combine  = [Population,EP];
    Selected = false(1,length(Combine));
    Selected(1:length(Population)) = true;
    Dis = pdist2(Combine.objs,Combine.objs);
    Dis(logical(eye(length(Dis)))) = inf;
    while sum(Selected) < min(N,length(Selected))
        subDis = sort(Dis(~Selected,Selected),2);
        max_value = max(prod(subDis(:,1:min(M,size(subDis,2))),2));
        max_indices = find(prod(subDis(:,1:min(M,size(subDis,2))),2) == max_value);
        random_index = max_indices(randi(length(max_indices)));
        Remain = find(~Selected);
        Selected(Remain(random_index)) = true;
    end
    % Add new subproblems
    newObjs = EP(Selected(length(Population)+1:end)).objs;
    W = [W;newObjs./repmat(sum(newObjs,2),1,size(newObjs,2))];
    % Add new solutions
    Population = Combine(Selected);
end