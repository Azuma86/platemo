isSaveCriticalDiff = 1;

% labels={'CAE-depth1','CAE-depth2','CAE'};
% % labels = {'CAE', 'CAEA', 'FTCA', 'TCA', 'SOINNplus', 'ASOINN'};
%labels = {'PPS','C3M','CAEAD','CMOEA-MS','CCMO','CM2T','CM2T-AWA'};
labels = {'ÉâÉìÉN','ç≈ëÂíl','å¬êî','çáåv'};
alpha=0.05; %?íòê´êÖïΩ0.1,0.05àΩ0.01
filename='RWMOP30~35.xlsx';
% % filename = 'statioanry_and_Nonstatioary_Comparison_AllAlgorithms_12datas.xlsx';
T=readtable(filename);
A = table2array(T);
%Acc_NMI_ARI_macro-F1=double(T);
%è∏èá
fig = criticaldifference_new(A,labels,alpha,0); 
%ç~èá
%fig = criticaldifference_new(1-A,labels,alpha,0); 

if isSaveCriticalDiff == 1
        filename = strcat('cd_nonstationary');
        SaveFig(fig, filename);
end




function SaveFig(ff, data_name)

drawnow;

temp.figunit = ff.Units;
ff.Units = 'centimeters';
ff.PaperPositionMode = 'manual';

% Change positions and paper size for removing wasted spaces in PDF
set(ff, 'PaperPosition', [0.3 -5.7 19.0 13.6]);
ff.PaperSize = [20.3, 6.6];

temp.figpaperunit = ff.PaperUnits;
ff.PaperUnits = 'centimeters';
temp.figsize = ff.PaperSize;

% print(ff,data_name,'-dpdf','-r300','-bestfit')
print(ff,'-vector','-r450',data_name,'-dpdf')
% print(ff,data_name,'-dpng','-r300')

ff.PaperSize = temp.figsize;
ff.Units = temp.figunit;
ff.PaperUnits = temp.figpaperunit;

end


