function [SubPopulation,SubPop_Constraint_value,SubPop_Scalarizing_value] = EnSele(Population,Pop_Constraint_value,Pop_Scalarizing_value,N)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort([Pop_Constraint_value Pop_Scalarizing_value],N);
    Next = FrontNo < MaxFNo;
    disp(MaxFNo)
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance([Pop_Constraint_value Pop_Scalarizing_value],FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    SubPopulation = Population(Next);
    SubPop_Constraint_value = Pop_Constraint_value(Next);
    SubPop_Scalarizing_value = Pop_Scalarizing_value(Next);
end