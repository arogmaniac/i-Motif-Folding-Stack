function FWDREVCD(filepathdown,filepathup,sampleID,CIToggle)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dirnameup=filepathdown;
dirnamedown=filepathup;

datafileupread=dir(fullfile(dirnameup,'*NumericalScatter.txt'));
datafileup=struct2cell(datafileupread);
datafiledownread=dir(fullfile(dirnamedown,'*NumericalScatter.txt'));
datafiledown=struct2cell(datafiledownread);

%%%%%%%%%%%%%%%%%%%%%Decreasing pH File Data Retrieval%%%%%%%%%%%%%%%%%%%%%
%Shows up as red trendline
%Read in the Data for the scatters
dfusp1=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [2,0,2,0]);
datausp1=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [3,0,2+dfusp1,1]);
nextrownum=2+dfusp1+2;

dfusp2=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [nextrownum,0,nextrownum,0]);
datausp2=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [nextrownum+1,0,nextrownum+dfusp2,1]);
nextrownum2=nextrownum+dfusp2+2;

dfusp3=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [nextrownum2,0,nextrownum2,0]);
datausp3=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [nextrownum2+1,0,nextrownum2+dfusp3,1]);
nextrownum3=nextrownum2+dfusp3+2;

%Read in the Data for the regression line
dfurl=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [nextrownum3,0,nextrownum3,0]);
dataurl=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [nextrownum3+1,0,nextrownum3+dfurl,1]);
nextrownum4=nextrownum3+dfurl+2;

%Read in the Data for the Upper and Lower CI
dfuCIu=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [nextrownum4,0,nextrownum4,0]);
datauCIu=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [nextrownum4+1,0,nextrownum4+dfuCIu,1]);
nextrownum5=nextrownum4+dfuCIu+2;
dfuCIl=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [nextrownum5,0,nextrownum5,0]);
datauCIl=dlmread(fullfile(filepathdown,datafileup{1,1}),'\t',...
    [nextrownum5+1,0,nextrownum5+dfuCIl,1]);

%%%%%%%%%%%%%%%%%%%%%Increasing pH File Data Retrieval%%%%%%%%%%%%%%%%%%%%%
%Shows up as black trend line and data points
%Read in the Data for the scatters
dflsp1=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [2,0,2,0]);
datalsp1=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [3,0,2+dflsp1,1]);
nextrownum=2+dflsp1+2;

dflsp2=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [nextrownum,0,nextrownum,0]);
datalsp2=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [nextrownum+1,0,nextrownum+dflsp2,1]);
nextrownum2=nextrownum+dflsp2+2;

dflsp3=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [nextrownum2,0,nextrownum2,0]);
datalsp3=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [nextrownum2+1,0,nextrownum2+dflsp3,1]);
nextrownum3=nextrownum2+dflsp3+2;

%Read in the Data for the regression line
dflrl=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [nextrownum3,0,nextrownum3,0]);
datalrl=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [nextrownum3+1,0,nextrownum3+dflrl,1]);
nextrownum4=nextrownum3+dflrl+2;

%Read in the Data for the Upper and Lower CI
dflCIu=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [nextrownum4,0,nextrownum4,0]);
datalCIu=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [nextrownum4+1,0,nextrownum4+dflCIu,1]);
nextrownum5=nextrownum4+dflCIu+2;
dflCIl=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [nextrownum5,0,nextrownum5,0]);
datalCIl=dlmread(fullfile(filepathup,datafiledown{1,1}),'\t',...
    [nextrownum5+1,0,nextrownum5+dflCIl,1]);

%%%%%%%%%%%%%%%%Generate Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
box on
hold on
axis([5 8 -0.1 1.1])
hold on
plot(dataurl(:,1),dataurl(:,2),'-r','LineWidth',2)
hold on
plot(datalrl(:,1),datalrl(:,2),'-k','LineWidth',2)
hold on
if CIToggle==1
    plot(datauCIu(:,1),datauCIu(:,2),'--r','LineWidth',2)
    hold on
    plot(datauCIl(:,1),datauCIl(:,2),'--r','LineWidth',2)
    hold on
    plot(datalCIu(:,1),datalCIu(:,2),'--k','LineWidth',2)
    hold on
    plot(datalCIl(:,1),datalCIl(:,2),'--k','LineWidth',2)
    hold on
else if CIToggle==0
    end
end
scatter(datausp1(:,1),datausp1(:,2),'or','filled','MarkerEdgeColor','k','Linewidth',1)
hold on
scatter(datausp2(:,1),datausp2(:,2),'or','filled','MarkerEdgeColor','k','Linewidth',1)
hold on
scatter(datausp3(:,1),datausp3(:,2),'or','filled','MarkerEdgeColor','k','Linewidth',1)
hold on
scatter(datalsp1(:,1),datalsp1(:,2),'ok','filled','MarkerEdgeColor','k','Linewidth',1)
hold on
scatter(datalsp2(:,1),datalsp2(:,2),'ok','filled','MarkerEdgeColor','k','Linewidth',1)
hold on
scatter(datalsp3(:,1),datalsp3(:,2),'ok','filled','MarkerEdgeColor','k','Linewidth',1)
hold on
set(get(gca,'XLabel'),'String','pH','FontName',...
    'Times New Roman','FontSize',36,'FontWeight','bold')
set(get(gca,'YLabel'),'String',...
    '{\itx}_{iM}','FontName',...
    'Times New Roman','FontSize',36,'Fontweight','bold')
set(gca,'XTick',5:.5:8)
set(gca,'XTickLabel',char('5.0','5.5','6.0','6.5',...
    '7.0','7.5','8.0'),'FontSize',28,'FontName','Times New Roman')
set(gca,'YTick',0:0.2:1)
set(gca,'YTickLabel',char('0.0','0.2','0.4','0.6','0.8','1.0'),...
    'FontSize',28,'FontName','Times New Roman')

graphname=strcat('C:\Users\Robert\Desktop','/',sampleID,'_FWDREVLowRes');
print(fullfile(graphname),'-dtiff','-r200');

graphname1=strcat('C:\Users\Robert\Desktop','/',sampleID,'_FWDREVHighRes');
print(fullfile(graphname1),'-dtiff','-r1000');














end