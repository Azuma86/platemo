classdef CTSEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained>
% Constrained two-stage evolutionary algorithm
%evaluation --- 1 ---final population of EP

%------------------------------- Reference --------------------------------
% F. Ming, W. Gong, H. Zhen, S. Li, L. Wang, and Z. Liao, A simple
% two-stage evolutionary algorithm for constrained multi-objective
% optimization, Knowledge-Based Systems, 2021, 228: 107263.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Fei Ming

    methods
        function main(Algorithm,Problem)
            %% Generate the sampling points and random population
            evaluation = Algorithm.ParameterSet(1);
            Population = Problem.Initialization();
            W = UniformPoint(Problem.N,Problem.M);
            [ARMOEA_Archive,RefPoint,Range] = UpdateRefPoint(Population.objs,W,[]);
            CV = sum(max(0,Population.cons),2);
            Archive = Population(CV==0);
            stage_conver = 0;
            EP = updateEP2(Population);
            
            %% Optimization
            while Algorithm.NotTerminated(Population)
                if Problem.FE<0.5*Problem.maxFE
                    % evolve population to PF by ARMOEA
                    MatingPool = MatingSelection1(Population,RefPoint,Range);
                    Offspring = OperatorGA(Problem,Population(MatingPool));
                    [ARMOEA_Archive,RefPoint,Range] = UpdateRefPoint([ARMOEA_Archive;Offspring.objs],W,Range);
                    Archive = UpdateArchive(Archive,[Population,Offspring],Problem.N);
                    [Population,Range] = EnvironmentalSelection1([Population,Offspring],RefPoint,Range,Problem.N);
                else
                    if stage_conver==0
                        % exchange archive and population
                        temp = Population;
                        Population = Archive;
                        Archive = temp;
                        stage_conver = 1;
                    end
                    % evolve population to CPF by modified SPEA2
                    MatingPool = MatingSelection2(Population,Archive,Problem.N);
                    Offspring = OperatorGA(Problem,MatingPool);
                    Population = EnvironmentalSelection2([Population,Offspring],Problem.N);
                end
                EP = updateEP2([EP,Offspring]);
                % Output the non-dominated and feasible solutions.
                if Problem.FE >= Problem.maxFE
                    if evaluation == 2
                        Population = EP;
                    end
                end
            end
        end
    end
end