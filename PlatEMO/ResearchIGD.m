function ResearchIGD()
    IGDList = [];
    runNo = 21;
    for i = 1:runNo
        load(['CM2TAWA2_LIRCMOP1_M2_D30_',num2str(i),'.mat']);
        IGD = metric.IGD(1);%Number of result
        IGDList = [IGDList; i,IGD];
    end
    sortedIGDList = sortrows(IGDList, 2);
    %disp(sortedIGDList(11,:));%中央値の番号指定
    disp(sortedIGDList(11,1));
    t=sortedIGDList(11,1);
    writematrix(t,'EP.txt','WriteMode','append')
end

