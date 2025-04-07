function [SubPopulation,SubPop_Constraint_value,SubPop_Scalarizing_value] = EpsilonEnvironmentalSelectionAWA(Population,Pop_Constraint_value,Pop_Scalarizing_value,N,epsilon)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% con check
    feasible_indices = find(Pop_Constraint_value<=epsilon);
    infeasible_indices = find(Pop_Constraint_value>epsilon);

    if length(feasible_indices)>=N
        Pop_Constraint_value(infeasible_indices)=Inf;
        Pop_Scalarizing_value(infeasible_indices)=Inf;
        % Non-dominated sorting
        [FrontNo,MaxFNo] = NDSort([Pop_Constraint_value Pop_Scalarizing_value],N);
        Next = FrontNo < MaxFNo;
        
        % Calculate the crowding distance of each solution
        CrowdDis = CrowdingDistance([Pop_Constraint_value Pop_Scalarizing_value],FrontNo);
        
        % Select the solutions in the last front based on their crowding distances
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(CrowdDis(Last),'descend');
        Next(Last(Rank(1:N-sum(Next)))) = true;
    else
        %% Non-dominated sorting
        Mat_Con_and_Scalar = [Pop_Constraint_value,Pop_Scalarizing_value];
        % 行番号を取得
        row_numbers = (1:size(Mat_Con_and_Scalar, 1))';
        % 行列と行番号を結合
        Mat_Con_and_Scalar_with_row_numbers = [Mat_Con_and_Scalar, row_numbers];
        SortedMat_Con_and_Scalar = sortrows(sortrows(Mat_Con_and_Scalar_with_row_numbers,2),1);
        Next = SortedMat_Con_and_Scalar(1:N,3);
    end
    
    %% Population for next generation
    SubPopulation = Population(Next);
    SubPop_Constraint_value = Pop_Constraint_value(Next);
    SubPop_Scalarizing_value = Pop_Scalarizing_value(Next);
    %FrontNo    = FrontNo(Next);
    %CrowdDis   = CrowdDis(Next);
end