function Plot_EP()
    IGDList = [];
    runNo = 21;
    numResult = 1;
    proNo = 9;
    M=2;
    D=15;
    algo = "CMOEA_MS_EP";
    pro = "MW";
    for i = 1:runNo
        load(sprintf('%s_%s%d_M%d_D%d_%d.mat',algo,pro,proNo,M,D,i));
        IGD = metric.IGD(numResult);%Number of result
        IGDList = [IGDList; i,IGD];
    end
    sortedIGDList = sortrows(IGDList, 2);
    ind = sortedIGDList(11,1);
    disp(sortedIGDList(11,1));%中央値の番号指定
    disp(sortedIGDList(11,2));
    load(sprintf('%s_%s%d_EP_M%d_D%d_%d.mat',algo,pro,proNo,M,D,ind),'EP');
    Obj = EP.objs;
    if M == 2
        scatter(Obj(:,1),Obj(:,2),40,[0.8 0.8 1],'filled','LineWidth',1.5,'MarkerEdgeColor','red','MarkerFaceColor','red','MarkerFaceAlpha','0.1');
    else
        scatter3(Obj(:,1),Obj(:,2),Obj(:,3),25,[0.8 0.8 1],'filled','LineWidth',1,'MarkerEdgeColor','red','MarkerFaceColor','red','MarkerFaceAlpha','0.05');
    end
end

