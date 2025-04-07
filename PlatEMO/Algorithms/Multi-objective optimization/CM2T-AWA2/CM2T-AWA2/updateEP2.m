function EP = updateEP2(EP)
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
    feasible_index = find(sum(max(0,EP.cons),2) == 0);
    EP = EP(feasible_index);
    EP = EP(NDSort(EP.objs,1)==1);
end