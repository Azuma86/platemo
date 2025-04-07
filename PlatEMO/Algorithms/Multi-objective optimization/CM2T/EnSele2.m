function [SubPopulation,SubPop_Constraint_value,SubPop_Scalarizing_value] = EnvironmentalSelection2(Population,Pop_Constraint_value,Pop_Scalarizing_value,N)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    [~,Rank] = sort(Pop_Scalarizing_value,'ascend');
    
    %% Population for next generation
    SubPopulation = Population(Rank(1:N));
    SubPop_Constraint_value = Pop_Constraint_value(Rank(1:N));
    SubPop_Scalarizing_value = Pop_Scalarizing_value(Rank(1:N));
end