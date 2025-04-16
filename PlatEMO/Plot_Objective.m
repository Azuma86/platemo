function Plot_Objective()
    IGDList = [];
    runNo = 21;
    numResult = 1;
    for i = 1:runNo
        load(['CM2T_RWMOP1_M2_D4_',num2str(i),'.mat']);
        IGD = metric.IGD(numResult);%Number of result
        IGDList = [IGDList; i,IGD];
    end
    sortedIGDList = sortrows(IGDList, 2);
    ind = sortedIGDList(11,1);
    disp(sortedIGDList(11,1));%中央値の番号指定
    load(sprintf('CM2T_RWMOP1_M2_D4_%d.mat',ind),'result');
    Population = result{numResult,2};
    disp(size(Population.objs))
    disp(size(Population.decs))
    Obj = result{numResult,2}.objs;
    feasible_index = find(sum(max(0,Population.cons),2) == 0);
    % 特定の行以外を抽出する
    infeasible_index = setdiff(1:size(Obj, 1), feasible_index);
    
    scatter(Obj(infeasible_index,1),Obj(infeasible_index,2),25,[0.8 0.8 1],'filled','MarkerEdgeColor','blue');
    hold on
    scatter(Obj(feasible_index,1),Obj(feasible_index,2),25,[1 0.8 0.8],'filled','MarkerEdgeColor','red');
    hold off
end

