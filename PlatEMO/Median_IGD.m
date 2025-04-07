function Median_IGD()
    runNo = 21;
    numResult = 1;
    proNo = 10;
    M=2;
    D=15;
    algo = ["PPS_EP","CMOEA_MS_EP","CCMO_EP","C3M_EP","CM2T_EP","CM2TAWA2"];
    pro = "MW";
    for j = 1:6
        IGDList = [];
        IGD = [];
        for i = 1:runNo
            load(sprintf('%s_%s%d_M%d_D%d_%d.mat',algo(j),pro,proNo,M,D,i));
            IGD = metric.IGD(numResult);%Number of result
            IGDList = [IGDList; i,IGD];
        end
        sortedIGDList = sortrows(IGDList, 2);
        ind = sortedIGDList(11,1);
        writematrix(sortedIGDList(10,2),'IGD.txt','WriteMode','append');%中央値の番号指定
    end
end