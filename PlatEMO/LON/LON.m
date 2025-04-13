pro = "RWMOP6";
m = 2;
d = 7;
runNo = 31;
k        = 10;      % k近傍
numStart = 1000;     % ランダムに選ぶ開始点数 
%NSGAIICVM,'NSGAIICDP','NSGAIINSGAIIchev','NSGAIINSGAIIcount'
algo = {'NSGAIICDP'};  
solutions = all_solution(algo,pro,m,d,runNo);
neighbor_search(solutions,pro,k,numStart);
