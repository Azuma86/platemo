pro = "RWMOP22";
m = 2;
d = 9;
runNo = 31;
%NSGAIICVM,'NSGAIINSGAIICDP','NSGAIINSGAIIchev','NSGAIINSGAIIcount'
algo = {'NSGAIICVM'};  
solutions = all_solution(algo,pro,m,d,runNo);
neighbor_search(solutions,pro);
