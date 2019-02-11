function PromptMasterCD(filepath,epsilonLumolpercm,pathlengthcm,...
    sampleID)
%This function will output several important data processing results from
%CD spectra pH Titrations of the DNA i-Motif or the G-quadruplex (if salt
%is the titrant). Essentially, triplicate data named
%"pH*100.replicatenumber.txt" are compiled, blanked using a file
%"Blank.txt", and then normalized. Following normalization, each titration
%data set is plotted on the same graph and a regression is performed to fit
%each set to a sigmoidal titration curve. The result is an output of a
%figure with each curve and data set in the same graph for visual
%comparison as well as a set of results to describe the parameters for the
%fitting routine including R^2, the fit values of which one is the pH
%transition value (defined as the pH at which half the population of DNA in
%the sample resides in the folded state), the error of each parameter and
%the concentrations of each replicate sample and their respective standard
%errors. The second part will then use a cubic spline fitting routine to 
%model each replicate data set and then plot those derivatives resulting in
%a graph made up of gaussian type bell curves. Since this data typically
%results in two distinct transitions, a second nonlinear regression is used
%to model the data in this derivative plot using the sum of two gaussian
%functions. As was performed for the sigmoidal plot, the parameters of the
%fit for each replicate set will be output along with their errors and R^2
%values in the results file. The code will compute the error propagation
%for reporting a single value for each transition as well as a compilation
%of the CD spectra for one of the replicate titrations to determine the
%wavelength for analysis. For the gaussian regression plot a plot of the
%residuals will be output as well to determine if the fit is appropriate
%and indeed models the data well. Finally, a plot of the UV difference
%spectra will also be output and the results for that will be the maximum
%and minimum UV wavelengths. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Section 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import the raw data from the file in which it is contained, convert into
%matrix format and assign indices for use later on in the data processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirName=filepath;            
files1=dir(fullfile(dirName,'*.1.txt'));
files2=dir(fullfile(dirName,'*.2.txt'));
files3=dir(fullfile(dirName,'*.3.txt'));
Blank=dir(fullfile(dirName,'Blank.txt'));

filenames1=struct2cell(files1);
filenames2=struct2cell(files2);
filenames3=struct2cell(files3);
Blankname=struct2cell(Blank);

file1=transpose(filenames1(1,:));
file2=transpose(filenames2(1,:));
file3=transpose(filenames3(1,:));
Blankready=transpose(Blankname(1,:));

n1=length(file1);
n2=length(file2);
n3=length(file3);

%Read in the cd and uv data from each file
file1cddata=zeros(1201,n1);
file1uvdata=zeros(1201,n1);
for i=1:n1;
    file1cddata(:,i)=dlmread(fullfile(filepath,file1{i,1}),'\t',...
        [20,1,1220,1]);
    file1uvdata(:,i)=dlmread(fullfile(filepath,file1{i,1}),'\t',...
        [20,2,1220,2]);
end

file2cddata=zeros(1201,n2);
file2uvdata=zeros(1201,n2);
for i=1:n2;
    file2cddata(:,i)=dlmread(fullfile(filepath,file2{i,1}),'\t',...
        [20,1,1220,1]);
    file2uvdata(:,i)=dlmread(fullfile(filepath,file2{i,1}),'\t',...
        [20,2,1220,2]);
end

file3cddata=zeros(1201,n3);
file3uvdata=zeros(1201,n3);
for i=1:n3;
    file3cddata(:,i)=dlmread(fullfile(filepath,file3{i,1}),'\t',...
        [20,1,1220,1]);
    file3uvdata(:,i)=dlmread(fullfile(filepath,file3{i,1}),'\t',...
        [20,2,1220,2]);
end

Blankcddata(:,1)=dlmread(fullfile(filepath,Blankready{1,:}),'\t',...
    [20,1,1220,1]);
Blankuvdata(:,1)=dlmread(fullfile(filepath,Blankready{1,:}),'\t',...
    [20,2,1220,2]);

%Sizes of matrices with data for index assignments Note that cd and uv data
%are the same size so no need to have different indexing for each set

[q1,r1]=size(file1cddata);
[q2,r2]=size(file2cddata);
[q3,r3]=size(file3cddata);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Section 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Blank the Data in both the UV and CD spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file1cddatacorr=zeros(q1,r1);
file1uvdatacorr=zeros(q1,r1);
for i=1:q1
    for j=1:r1
        file1cddatacorr(i,j)=file1cddata(i,j)-Blankcddata(i,1);
        file1uvdatacorr(i,j)=file1uvdata(i,j)-Blankuvdata(i,1);
    end
end

file2cddatacorr=zeros(q2,r2);
file2uvdatacorr=zeros(q2,r2);
for i=1:q2
    for j=1:r2
        file2cddatacorr(i,j)=file2cddata(i,j)-Blankcddata(i,1);
        file2uvdatacorr(i,j)=file2uvdata(i,j)-Blankuvdata(i,1);
    end
end

file3cddatacorr=zeros(q3,r3);
file3uvdatacorr=zeros(q3,r3);
for i=1:q3
    for j=1:r3
        file3cddatacorr(i,j)=file3cddata(i,j)-Blankcddata(i,1);
        file3uvdatacorr(i,j)=file3uvdata(i,j)-Blankuvdata(i,1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%Section 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalize the CD spectra for the transition plot based on the absorbance
%at 260 nm denoted as row 601 in the imported UV spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Extract the absorbance at 260nm from each spectra
file1uv260=zeros(1,r1);
for i=1:r1;
    file1uv260(1,i)=file1uvdatacorr(601,i);
end
file1uvavg=mean(file1uv260);
file1uvstd=std(file1uv260);

file2uv260=zeros(1,r2);
for i=1:r2;
    file2uv260(1,i)=file2uvdatacorr(601,i);
end
file2uvavg=mean(file2uv260);
file2uvstd=std(file2uv260);

file3uv260=zeros(1,r3);
for i=1:r3;
    file3uv260(1,i)=file3uvdatacorr(601,i);
end
file3uvavg=mean(file3uv260);
file3uvstd=std(file3uv260);

%Determine the concentration and std dev using the epislon value in 
%(L umol cm^-1)

file1conc=file1uvavg/(epsilonLumolpercm*pathlengthcm);
file1concstd=file1uvstd/(epsilonLumolpercm*pathlengthcm);

file2conc=file2uvavg/(epsilonLumolpercm*pathlengthcm);
file2concstd=file2uvstd/(epsilonLumolpercm*pathlengthcm);

file3conc=file3uvavg/(epsilonLumolpercm*pathlengthcm);
file3concstd=file3uvstd/(epsilonLumolpercm*pathlengthcm);

concdata=[file1conc,file2conc,file3conc;...
    file1concstd,file2concstd,file3concstd];

%Normalize the cd spectra at the wavelength maximum indicated by "row" 
%using the UV260nm concentration average for all spectra to plot the 
%titration curves

[~,row]=max(file1cddatacorr(1:1200,1));

file1cdtrans=zeros(1,r1);
for i=1:r1;
    file1cdtrans(1,i)=file1cddatacorr(row,i)./(1000000*file1conc*pathlengthcm);
end
max1=max(file1cdtrans);
min1=min(file1cdtrans);
file1cdnorm=zeros(1,r1);
for i=1:r1;
    file1cdnorm(1,i)=(file1cdtrans(1,i)-min(file1cdtrans))/...
        (max(file1cdtrans)-min(file1cdtrans));
end

file2cdtrans=zeros(1,r2);
for i=1:r2;
    file2cdtrans(1,i)=file2cddatacorr(row,i)./(1000000*file2conc*pathlengthcm);
end
max2=max(file2cdtrans);
min2=min(file2cdtrans);
file2cdnorm=zeros(1,r2);
for i=1:r2;
    file2cdnorm(1,i)=(file2cdtrans(1,i)-min(file2cdtrans))/...
        (max(file2cdtrans)-min(file2cdtrans));
end

file3cdtrans=zeros(1,r3);
for i=1:r3;
    file3cdtrans(1,i)=file3cddatacorr(row,i)./(1000000*file3conc*pathlengthcm);
end
max3=max(file3cdtrans);
min3=min(file3cdtrans);
file3cdnorm=zeros(1,r3);
for i=1:r3;
    file3cdnorm(1,i)=(file3cdtrans(1,i)-min(file3cdtrans))/...
        (max(file3cdtrans)-min(file3cdtrans));
end

%Create the x values (pH) from the filenames which give the pH values
phvals1(1,:)=cell2mat(filenames1(1,:));
phvals2(1,:)=cell2mat(filenames2(1,:));
phvals3(1,:)=cell2mat(filenames3(1,:));

phvals1cell{1,r1}=[];
for i=1:r1-1;
    phvals1cell{1,1}=phvals1(1,1:3);
    phvals1cell{1,i+1}=phvals1(1,(i*9+1):(i*9+1)+2);
end

phvals2cell{1,r2}=[];
for i=1:r2-1;
    phvals2cell{1,1}=phvals2(1,1:3);
    phvals2cell{1,i+1}=phvals2(1,(i*9+1):(i*9+1)+2);
end

phvals3cell{1,r3}=[];
for i=1:r3-1;
    phvals3cell{1,1}=phvals3(1,1:3);
    phvals3cell{1,i+1}=phvals3(1,(i*9+1):(i*9+1)+2);
end

phvals1mat=str2double(phvals1cell);
phvals2mat=str2double(phvals2cell);
phvals3mat=str2double(phvals3cell);

%Create an X,Y set for the transitions and concatenate the vectors 

phtransitionplot1=[transpose(phvals1mat).*0.01,transpose(file1cdnorm)];
phtransitionplot2=[transpose(phvals2mat).*0.01,transpose(file2cdnorm)];
phtransitionplot3=[transpose(phvals3mat).*0.01,transpose(file3cdnorm)];
[phx1,~]=size(phtransitionplot1);
[phx2,~]=size(phtransitionplot2);
[phx3,~]=size(phtransitionplot3);

phtransitionplotALL=vertcat(phtransitionplot1,phtransitionplot2,...
    phtransitionplot3);

%%%%%%%%%%%%%%%%%%%%%%%%%Section 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sigmoidal Regression (global) of the normalized replicate titration data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Perform the regression analysis for each set of data

modelfun=@(b,x)((-1./(1+exp((-b(1).*x)+b(2))))+1);
betazero=[5;32.5];
[beta1,R1,~,CovB1]=nlinfit(phtransitionplotALL(:,1),...
    phtransitionplotALL(:,2),modelfun,betazero);

%Determine values to plot each regression line

regressionvalsx1=transpose(3.5:0.01:8.5);
[q]=length(regressionvalsx1);
regressionvalsy1=zeros(q,1);
for i=1:q;
    regressionvalsy1(i,1)=((-1./(1+exp((-beta1(1,1).*regressionvalsx1(i,1))...
        +beta1(2,1))))+1);
end

%Calculate the 95% CI values for each regression for graphical use
[Ypred,delta] = nlpredci(modelfun,regressionvalsx1,beta1,R1,'Covar',CovB1);

Ypredupper2=Ypred+delta;
Ypredlower2=Ypred-delta;

%Calculate the fit statistics.
transition1=beta1(2,1)/beta1(1,1);
transition1stddev=transition1*sqrt(((sqrt(CovB1(2,2))/beta1(2,1))^2)+...
    ((sqrt(CovB1(1,1))/beta1(1,1))^2));
beta1stddev=sqrt(CovB1(1,1));
beta2stddev=sqrt(CovB1(2,2));
meany=mean(phtransitionplotALL(:,2));

SStotal=(phtransitionplotALL(:,2)-meany)'*(phtransitionplotALL(:,2)-meany);
SSresid=R1'*R1;

Rsquared=1-(SSresid/SStotal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Section 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figures and statistics for the global sigmoidal regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot the transitions from each replicate experiment on one graph

figure(1)
box on
hold on
plot(regressionvalsx1(:,1),Ypredupper2(:,1),'--r','LineWidth',2,'Color','k');
hold on
legend({'\Delta_{95%}'},'FontName','Times New Roman','FontSize',24)
hold on
plot(regressionvalsx1(:,1),Ypredlower2(:,1),'--r','LineWidth',2,'Color','k')
hold on
plot(regressionvalsx1(:,1),regressionvalsy1(:,1),'-k','LineWidth',2,...
    'Color','k');
hold on
scatter(phtransitionplotALL(1:phx1,1),phtransitionplotALL(1:phx1,2),35,...
    'or','filled','MarkerEdgeColor','k','Linewidth',2)
hold on
scatter(phtransitionplotALL(phx1+1:(phx1+phx2),1),...
    phtransitionplotALL((phx1+1):(phx1+phx2),2),35,'or','filled',...
    'MarkerEdgeColor','k','Linewidth',2)
hold on
scatter(phtransitionplotALL((phx1+phx2+1):(phx1+phx2+phx3),1),...
    phtransitionplotALL((phx1+phx2+1):(phx1+phx2+phx3),2),35,'or',...
    'filled','MarkerEdgeColor','k','Linewidth',2)
axis([5 8 -0.1 1.1])
set(get(gca,'XLabel'),'String','pH','FontName',...
    'Times New Roman','FontSize',32,'FontWeight','bold')
set(get(gca,'YLabel'),'String',...
    '{\itx}_{iM}','FontName',...
    'Times New Roman','FontSize',32,'Fontweight','bold')
set(gca,'XTick',5:0.5:8)
set(gca,'XTickLabel',char('5.0','5.5','6.0','6.5',...
    '7.0','7.5','8.0'),'FontSize',24,'FontName','Times New Roman')
set(gca,'YTick',0:0.2:1)
set(gca,'YTickLabel',char('0.0','0.2','0.4','0.6','0.8','1.0'),...
    'FontSize',24,'FontName','Times New Roman')

graphname1=strcat(filepath,strcat('/',sampleID,'_SigmoidpHTPlotLowRes'));
print(fullfile(graphname1),'-dtiff','-r200');

graphname2=strcat(filepath,strcat('/',sampleID,'_SigmoidpHTPlotHighRes'));
print(fullfile(graphname2),'-dtiff','-r1000');

%%%Prompt the user to guess the modelfuntoggle (i.e. 1,2,or 3 gauss) and
%%%then write the initial guess values for the locations of each peak based
%%%on the raw data plotted in figure 1)

modelfuntoggle=...
    input('type the model function guess i.e. 1,2 or 3 Gauss\n');
initialguesstoggle=input(strcat(...
    'type the vector with init. guesses for b vals\n',...
    ' i.e. [height pH width]\n',...
    'always use 1 for height and go up in tenths from 0.1 for width\n',...
    '***for width, 0.1 means the transition occurs over 1 pH unit***\n'));
fprintf('\n WORKING NOW \n');

%Print a single compiled CD spectra with normalization of the spectra from
%###.1.txt data set

%Normalize whole CD spectra set
cdnormplot1=zeros(q1,r1);
for i=1:r1;
    for j=1:q1;
        cdnormplot1(j,i)=file1cddatacorr(j,i)./(1*file1conc*pathlengthcm);
    end
end

UVwavelength=dlmread(fullfile(filepath,file1{1,1}),'\t',[20,0,1220,0]);
Plotmat=[UVwavelength,cdnormplot1];

%This section creates data to give reference axes internally within the
%compiled CD spectra plot at theta=0 and at cdwavenm. 
horzgridlinesx=200:320;
horzgridlinesy=zeros(1,121);

%This section determines the coloration for the plots to cycle through. As
%long as there are not more than 20 pH values, you will not see a
%color repeat. The first file you input will be red which will gradually
%transition to blue which will gradually transition to green.
kolormap=zeros(r1,3);
midpt=ceil(r1/2);
for i=1:midpt;
    kolormap(i,1)=1-(i/(midpt));
    kolormap(i,2)=0;
    kolormap(i,3)=0+(i/(midpt));
end
for i=(midpt+1):r1;
    kolormap(i,1)=0;
    kolormap(i,2)=1-((r1-i)/i);
    kolormap(i,3)=0+((r1-i)/i);
end

%This section makes the CD spectra compiled figure for all pH values, salt
%concentration values, temperature values or whatever parameter you named
%your files for.

figure(2)
box on
hold on
plot(horzgridlinesx,horzgridlinesy,':k','LineWidth',2)
hold on
set(gca,'NextPlot','replacechildren','ColorOrder',kolormap)
hold on
plot(Plotmat(:,1),Plotmat(:,2:r1),'LineWidth',2);
legend off
set(gcf,'Colormap',flipud(kolormap(1:r1-1,:)))
axis([200 320 -8 16])
set(get(gca,'XLabel'),'String','\lambda (nm)','FontName',...
    'Times New Roman','FontSize',32,'FontWeight','bold')
set(get(gca,'YLabel'),'String',...
    '[\theta] x 10^{6}(deg cm^{2} dmol^{-1})','FontName',...
    'Times New Roman','FontSize',32,'Fontweight','bold')
set(gca,'XTick',200:40:320)
set(gca,'XTickLabel',{'200','240','280','320'},...
    'FontSize',24,'FontName','Times New Roman')
set(gca,'YTick',-8:4:16)
set(gca,'YTickLabel',{' ',' -4','  0','  4','  8',' 12',' 16'},...
    'FontSize',24,'FontName','Times New Roman')
hcb=colorbar('location','EastOutside','YTick',...
    1:length(phtransitionplot1),'YTickLabel',...
    sprintf('%1.2f|',flipud(phtransitionplot1(:,1))),...
    'FontName','Times New Roman','FontSize',18);
set(get(hcb,'Title'),'String','pH','FontName',...
    'Times New Roman','FontSize',24,'FontWeight','bold')

graphname3=strcat(filepath,strcat('/',sampleID,...
    '_CDSpectraCompilationLowRes'));
print(fullfile(graphname3),'-dtiff','-r200');

graphname4=strcat(filepath,strcat('/',sampleID,...
    '_CDSpectraCompilationHighRes'));
print(fullfile(graphname4),'-dtiff','-r1000');

fprintf('\n CD SPECTRA COMPLETE \n');

%Residuals plot
figure(3)
box on
hold on
plot(4.5:.01:8.5,zeros(size(4.5:0.01:8.5)),':k','LineWidth',2)
hold on
scatter(phtransitionplotALL(:,1),R1,35,'ow','filled',...
    'MarkerEdgeColor','k','LineWidth',2)
hold on
axis([5 8 -0.35 0.35])
set(get(gca,'XLabel'),'String','pH','FontName',...
    'Times New Roman','FontSize',32,'FontWeight','bold')
set(get(gca,'YLabel'),'String',...
    'Residuals','FontName',...
    'Times New Roman','FontSize',32,'Fontweight','bold')
set(gca,'XTick',5:.5:8)
set(gca,'XTickLabel',char('5.0','5.5','6.0','6.5',...
    '7.0','7.5','8.0'),'FontSize',24,'FontName','Times New Roman')
set(gca,'YTick',-0.3:0.1:0.3)
set(gca,'YTickLabel',char('-0.3','-0.2','-0.1',' 0.0',' 0.1',' 0.2',' 0.3'),...
    'FontSize',24,'FontName','Times New Roman')

graphname5=strcat(filepath,strcat('/',sampleID,'_SigmoidResidPlotLowRes'));
print(fullfile(graphname5),'-dtiff','-r200');

fprintf('\n SIGMOID RESIDUALS COMPLETE \n');

%print pertinent info from regressions to file
resultspath=strcat(filepath,strcat('/',sampleID,'_SigmoidResults.txt'));
fid=fopen(resultspath,'wt');

printsampleID=sampleID;
fprintf(fid,printsampleID);
printcd='CD Analysis Wavelength = %3.1f nm\n';
fprintf(fid,printcd,UVwavelength(row,1));
printuv='UV Analysis Wavelength = 260 nm\n';
fprintf(fid,printuv);
printepsilon='Molar Absorptivity = %.4f (L cm uM^{-1})\n';
fprintf(fid,printepsilon,epsilonLumolpercm);
printconcentrationdata='\n\n ----Concentration Values and Std. Dev.----';
fprintf(fid,printconcentrationdata);
printconcentrationwave='\n                      @ 260 nm               \n';
fprintf(fid,printconcentrationwave);
fprintf(fid,'concentration = %2.4f\t\t\t%1.4f\n',concdata);
prntregstat='\n\n ----Sigmoidal Regression Statistics----';
prntregstat1='\n\n Equation\n              -1\n';
prntregstat2='y= ----------------------- + 1\n   1+exp((-b(1).*x)+b(2))\n';
fprintf(fid,prntregstat);
fprintf(fid,prntregstat1);
fprintf(fid,prntregstat2);
printbeta='\nb(1) = %3.4f;  b(2) = %3.4f;';
fprintf(fid,printbeta,beta1);
printbetastddev='\nStd. Devs \nb(1) = +/- %3.4f; b(2) = +/- %3.4f;';
fprintf(fid,printbetastddev,[beta1stddev,beta2stddev]);

printtransition='\nThe Transition Value is %3.4f';
printtransitionstddev='\nStandard Deviation is = +/- %3.4f'; 
fprintf(fid,printtransition,transition1);
fprintf(fid,printtransitionstddev,transition1stddev);
printRsquared='\n R^2 = %3.4f';
fprintf(fid,printRsquared,Rsquared);

fclose(fid); %Closes results file 
fprintf('\n SIGMOID RESULTS FILE COMPLETE \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%Section 6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gaussian bell curve fitting for the replicate data based on the slope
%between each successive point in the replicate titration experiments. In
%order to provide a better fit, a cubic spline interpolation is used on the
%original sigmoidal data and then it's derivative is plotted which prevents
%erroneous fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phtnorm1=[transpose(phvals1mat).*0.01,transpose(file1cdnorm)];
phtnorm2=[transpose(phvals2mat).*0.01,transpose(file2cdnorm)];
phtnorm3=[transpose(phvals3mat).*0.01,transpose(file3cdnorm)];

i1=length(phtnorm1);
i2=length(phtnorm2);
i3=length(phtnorm3);

slope1=zeros(i1-1,2);
for i=1:i1-1
    slope1(i,1)=(phtnorm1(i,1)+phtnorm1(i+1,1))/2;
    slope1(i,2)=-((phtnorm1(i+1,2)-phtnorm1(i,2))/...
        (phtnorm1(i+1,1)-phtnorm1(i,1)));
end

slope2=zeros(i2-1,2);
for i=1:i2-1
    slope2(i,1)=(phtnorm2(i,1)+phtnorm2(i+1,1))/2;
    slope2(i,2)=-((phtnorm2(i+1,2)-phtnorm2(i,2))/...
        (phtnorm2(i+1,1)-phtnorm2(i,1)));
end

slope3=zeros(i3-1,2);
for i=1:i3-1
    slope3(i,1)=(phtnorm3(i,1)+phtnorm3(i+1,1))/2;
    slope3(i,2)=-((phtnorm3(i+1,2)-phtnorm3(i,2))/...
        (phtnorm3(i+1,1)-phtnorm3(i,1)));
end

%Form one single matrix of data containing slope at avg pH points. This
%creates the X,Y data for the Gaussian Bell Curve regression. All slope
%values from each of the three replicates are included in the data matrix
%for regression as one data set.

gaussdata=vertcat(slope1,slope2,slope3);

%%%%%%%%%%%%%%%%%%%%Normalize the slopes to a max of 1%%%%%%%%%%%%%%%%%%%%%
gaussdatanorm=zeros(i1+i2+i3-3,1);
for i=1:i1+i2+i3-3;
    gaussdatanorm(i,1)=gaussdata(i,2)/max(gaussdata(:,2));
end

normdata=[gaussdata(:,1),gaussdatanorm];
normdatasorted=sortrows(normdata,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%Sum of Gaussian Fit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create scatter plot of the slope(i,i+1) vs. pHavg(i,i+1) data with.


%Non-linear fitting routine to estimate the b(#) parameters to be fit
%with beta0 indicating the initial guesses. These initial guesses are
%set for this particular data set, estimated by observations during
%execution of the experiment. b(1) and b(4) are the maxima of the gaussian
%peaks, so set one of them to 1 in the intial guess; b(2) and b(5) are the
%transition location in pH units so set those to 
fprintf('\n WORKING ON GAUSS REGRESSION \n');

if modelfuntoggle==1
    modelfun2=@(b,x)((b(1)*exp(-((x-b(2)).^2)/b(3))));
    beta02=initialguesstoggle;
else if modelfuntoggle==2
        modelfun2=@(b,x)((b(1)*exp(-((x-b(2)).^2)/b(3)))...
            +(b(4)*exp(-((x-b(5)).^2)/b(6))));
        beta02=initialguesstoggle;
else if modelfuntoggle==3
        modelfun2=@(b,x)((b(1)*exp(-((x-b(2)).^2)/b(3)))...
            +(b(4)*exp(-((x-b(5)).^2)/b(6)))...
            +(b(7)*exp(-((x-b(8)).^2)/b(9))));
        beta02=initialguesstoggle;
    end
    end
end

[beta2,R2,~,CovB2]=nlinfit(normdatasorted(:,1),normdatasorted(:,2)...
    ,modelfun2,beta02);

Testdata2=(4:0.05:8)';

fprintf('\n GAUSS REGRESSION COMPLETE \n');

%Calculate area of each bell curve
%first, build the equation
if modelfuntoggle==1;
    farea1=@(x)(beta2(1,1)*exp(-((x-beta2(1,2)).^2)/beta2(1,3)));
    area1=integral(farea1,4.5,8);
else if modelfuntoggle==2;
        farea12=@(x)(beta2(1,1)*exp(-((x-beta2(1,2)).^2)/beta2(1,3)));
        area12=integral(farea12,4.5,8);
        farea22=@(x)(beta2(1,4)*exp(-((x-beta2(1,5)).^2)/beta2(1,6)));
        area22=integral(farea22,4.5,8);
        percentarea12=(area12/(area12+area22))*100;
        percentarea22=(area22/(area12+area22))*100;
    else if modelfuntoggle==3;
            farea13=@(x)(beta2(1,1)*exp(-((x-beta2(1,2)).^2)/beta2(1,3)));
            area13=integral(farea13,4.5,8);
            farea23=@(x)(beta2(1,4)*exp(-((x-beta2(1,5)).^2)/beta2(1,6)));
            area23=integral(farea23,4.5,8);
            farea33=@(x)(beta2(1,7)*exp(-((x-beta2(1,8)).^2)/beta2(1,9)));
            area33=integral(farea33,4.5,8);
            percentarea13=(area13/(area13+area23+area33))*100;
            percentarea23=(area23/(area13+area23+area33))*100;
            percentarea33=(area33/(area13+area23+area33))*100;
        end
    end
end

fprintf('\n POPULATION STATS COMPLETE \n');
%Generate values for regression plot and 95% CI plots
[Ypred2,delta2]=nlpredci(modelfun2,Testdata2,beta2,R2,'Covar',CovB2);
Ypredupper2=Ypred2+delta2;
Ypredlower2=Ypred2-delta2;

%Create the Sum of Gaussians pH Transition plot
figure(4)
box on
hold on
h=plot(Testdata2,Ypredupper2,'--r','LineWidth',2);
hold on
legend(h,{'\Delta_{95%}'},'FontName','Times New Roman','FontSize',36)
hold on
plot(Testdata2,Ypred2,'k','LineWidth',3)
hold on
plot(Testdata2,Ypredlower2,'--r','LineWidth',2)
hold on
scatter(normdatasorted(:,1),normdatasorted(:,2),35,'or','filled',...
    'MarkerEdgeColor','k','LineWidth',2);
hold on
axis([5 8 -0.05 1.1])
hold on
set(get(gca,'XLabel'),'String','pH','FontName',...
    'Times New Roman','FontSize',48,'FontWeight','bold')
set(get(gca,'YLabel'),'String','\Delta {\it{x}}/\Delta pH',...
    'FontName','Times New Roman','FontSize',48,'FontWeight','bold')
set(gca,'XTick',5:0.5:8)
set(gca,'XTickLabel',char('5.0','5.5','6.0','6.5',...
    '7.0','7.5','8.0'),'FontSize',36,'FontName','Times New Roman')
set(gca,'YTick',0:0.2:1)
set(gca,'YTickLabel',char('0.0','0.2','0.4','0.6','0.8','1.0'),...
    'FontSize',36,'FontName','Times New Roman')

graphname6=strcat(filepath,'/',sampleID,'_',num2str(modelfuntoggle),...
    '_GausspHTPlotLowRes');
print(fullfile(graphname6),'-dtiff','-r200');

graphname7=strcat(filepath,'/',sampleID,'_',num2str(modelfuntoggle),...
    '_GausspHTPlotHighRes');
print(fullfile(graphname7),'-dtiff','-r1000');

%PRINT A GAUSS REGRESSION COMPARISON FILE FOR USE IN TRANSCOMP

resultspath=strcat(filepath,'/',sampleID,'_',num2str(modelfuntoggle),...
    '_GaussPeakScatter.txt');
fid=fopen(resultspath,'wt');

fprintf(fid,sampleID);
fprintf(fid,'\n%3.4f\t%3.4f\t%3.4f\t%3.4f',...
    [Testdata2,Ypred2,Ypredupper2,Ypredlower2].');

%Create a residuals plot for the two gaussian pH Transition plot
%Important parameters for the report. Beta std. dev.'s calculated using
    %the variance-covariance matrix. 

fprintf('\n GAUSS REGRESSION PLOT COMPLETE \n');   
%GAUSS1
if modelfuntoggle==1
   [Beta1error]=sqrt(CovB2(1,1));
   [pHTransitionerror1]=sqrt(CovB2(2,2));
   [Beta3error]=sqrt(CovB2(3,3));
else if modelfuntoggle==2
        [Beta1error]=sqrt(CovB2(1,1));
        [pHTransitionerror1]=sqrt(CovB2(2,2));
        [Beta3error]=sqrt(CovB2(3,3));
        [Beta4error]=sqrt(CovB2(4,4));
        [pHTransitionerror2]=sqrt(CovB2(5,5));
        [Beta6error]=sqrt(CovB2(6,6));
    else if modelfuntoggle==3
            [Beta1error]=sqrt(CovB2(1,1));
            [pHTransitionerror1]=sqrt(CovB2(2,2));
            [Beta3error]=sqrt(CovB2(3,3));
            [Beta4error]=sqrt(CovB2(4,4));
            [pHTransitionerror2]=sqrt(CovB2(5,5));
            [Beta6error]=sqrt(CovB2(6,6));
            [Beta7error]=sqrt(CovB2(7,7));
            [pHTransitionerror3]=sqrt(CovB2(8,8));
            [Beta9error]=sqrt(CovB2(9,9));
        end
    end
end

meany2=mean(normdatasorted(:,2));

SStotal2=(normdatasorted(:,2)-meany2)'*(normdatasorted(:,2)-meany2);
SSresid2=R2'*R2;

Rsquared2=1-(SSresid2/SStotal2);

%%%%%%plot the residuals for sum of Gauss fit%%%%%%%

figure(5)
box on
hold on
scatter(normdatasorted(:,1),R2,35,'ow','filled','MarkerEdgeColor','k',...
    'LineWidth',2);
hold on
plot(4:1:8,[0,0,0,0,0],':k','LineWidth',2);
axis([5 8 -0.35 0.35]);
set(get(gca,'XLabel'),'String','pH','FontName',...
    'Times New Roman','FontSize',32,'FontWeight','bold')
set(get(gca,'YLabel'),'String',...
    'Residuals','FontName',...
    'Times New Roman','FontSize',32,'Fontweight','bold')
set(gca,'XTick',5:0.5:8)
set(gca,'XTickLabel',char('5.0','5.5','6.0','6.5',...
    '7.0','7.5','8.0'),'FontSize',24,'FontName','Times New Roman')
set(gca,'YTick',-0.3:0.1:0.3)
set(gca,'YTickLabel',char('-0.3','-0.2','-0.1',' 0.0',' 0.1',...
    ' 0.2',' 0.3'),'FontSize',24,'FontName','Times New Roman')

graphname8=strcat(filepath,'/',sampleID,'_',num2str(modelfuntoggle),...
    '_GaussResidPlotLowRes');
print(fullfile(graphname8),'-dtiff','-r200');

fprintf('\n GAUSS RESIDUALS COMPLETE \n');

%Results from Two Gauss Regression
resultspath=strcat(filepath,'/',sampleID,'_',num2str(modelfuntoggle),...
    '_GaussResults.txt');
fid=fopen(resultspath,'wt');
printsampleID=sampleID;
fprintf(fid,printsampleID);
printcd='\n\nCD Analysis Wavelength = %3.1f nm\n';
fprintf(fid,printcd,UVwavelength(row,1));
printuv='UV Analysis Wavelength = 260 nm\n';
fprintf(fid,printuv);
printepsilon='Molar Absorptivity = %.4f (L cm uM^{-1})\n';
fprintf(fid,printepsilon,epsilonLumolpercm);
printconcentrationdata='\n\n ----Concentration Values and Std. Dev.----';
fprintf(fid,printconcentrationdata);
printconcentrationwave='\n                      @ 260 nm               \n';
fprintf(fid,printconcentrationwave);
fprintf(fid,'%2.4f\t\t\t%1.4f\n',concdata);
prntregstat='\n\n ----Bell Curve Regression Statistics----';
prntregstat1='\n\n Equation\n\n';
prntregstat2=...
    'Delta[x] = Sum from n=1:3 of b(n)*exp^{(-((pH-pHT(n))^{2})/b(n+2))}';
fprintf(fid,prntregstat);
fprintf(fid,prntregstat1);
fprintf(fid,prntregstat2);
if modelfuntoggle==1
    printbetaset1='\n b1 = %3.4f;\n b2 = pHT1 = %3.4f;\n b3 = %3.4f';
    fprintf(fid,printbetaset1,beta2(1,1:3));
    printbeta1error='\n sigma b1 = +/- %3.4f';
    printbeta2error='\n sigma pHT1 = +/- %3.4f';
    printbeta3error='\n sigma b3 = +/- %3.4f';
    printarea1='\n area pHT1 = %3.4f';
    fprintf(fid,printbeta1error,Beta1error);
    fprintf(fid,printbeta2error,pHTransitionerror1);
    fprintf(fid,printbeta3error,Beta3error);
    fprintf(fid,printarea1,area1);
else if modelfuntoggle==2
        printbetaset1='\n b1 = %3.4f;\n b2 = pHT1 = %3.4f;\n b3 = %3.4f';
        printbetaset2='\n b4 = %3.4f;\n b5 = pHT2 = %3.4f;\n b6 = %3.4f';
        fprintf(fid,printbetaset1,beta2(1,1:3));
        fprintf(fid,printbetaset2,beta2(1,4:6));
        printbeta1error='\n sigma b1 = +/- %3.4f';
        printbeta2error='\n sigma pHT1 = +/- %3.4f';
        printbeta3error='\n sigma b3 = +/- %3.4f';
        printbeta4error='\n sigma b4 = +/- %3.4f';
        printbeta5error='\n sigma pHT2 = +/- %3.4f';
        printbeta6error='\n sigma b6 = +/- %3.4f';
        printarea12='\n area percent pHT1 = %3.4f';
        printarea22='\n area percent pHT2 = %3.4f';
        fprintf(fid,printbeta1error,Beta1error);
        fprintf(fid,printbeta2error,pHTransitionerror1);
        fprintf(fid,printbeta3error,Beta3error);
        fprintf(fid,printbeta4error,Beta4error);
        fprintf(fid,printbeta5error,pHTransitionerror2);
        fprintf(fid,printbeta6error,Beta6error);
        fprintf(fid,printarea12,percentarea12);
        fprintf(fid,printarea22,percentarea22);
    else if modelfuntoggle==3
            printbetaset1='\n b1 = %3.4f;\n b2 = pHT1 = %3.4f;\n b3 = %3.4f';
            printbetaset2='\n b4 = %3.4f;\n b5 = pHT2 = %3.4f;\n b6 = %3.4f';
            printbetaset3='\n b7 = %3.4f;\n b8 = pHT3 = %3.4f;\n b9 = %3.4f';
            fprintf(fid,printbetaset1,beta2(1,1:3));
            fprintf(fid,printbetaset2,beta2(1,4:6));
            fprintf(fid,printbetaset3,beta2(1,7:9));
            printbeta1error='\n sigma b1 = +/- %3.4f';
            printbeta2error='\n sigma pHT1 = +/- %3.4f';
            printbeta3error='\n sigma b3 = +/- %3.4f';
            printbeta4error='\n sigma b4 = +/- %3.4f';
            printbeta5error='\n sigma pHT2 = +/- %3.4f';
            printbeta6error='\n sigma b6 = +/- %3.4f';
            printbeta7error='\n sigma b7 = +/- %3.4f';
            printbeta8error='\n sigma pHT3 = +/- %3.4f';
            printbeta9error='\n sigma b9 = +/- %3.4f';
            printarea13='\n area percent pHT1 = %3.4f';
            printarea23='\n area percent pHT2 = %3.4f';
            printarea33='\n area percent pHT3 = %3.4f';
            fprintf(fid,printbeta1error,Beta1error);
            fprintf(fid,printbeta2error,pHTransitionerror1);
            fprintf(fid,printbeta3error,Beta3error);
            fprintf(fid,printbeta4error,Beta4error);
            fprintf(fid,printbeta5error,pHTransitionerror2);
            fprintf(fid,printbeta6error,Beta6error);
            fprintf(fid,printbeta7error,Beta7error);
            fprintf(fid,printbeta8error,pHTransitionerror3);
            fprintf(fid,printbeta9error,Beta9error);
            fprintf(fid,printarea13,percentarea13);
            fprintf(fid,printarea23,percentarea23);
            fprintf(fid,printarea33,percentarea33);
        end
    end
end
printRsquared='\n\nR^2 = %3.4f';
fprintf(fid,printRsquared,Rsquared2);
fclose(fid);

fprintf('\n GUASS RESULTS FILE COMPLETE \n');

%%%%%%%%%%%%%%%%%%%%%Isothermal UV Difference Spectra (IDS)%%%%%%%%%%%%%%%%
dirName=filepath;            
pullfiles=dir(fullfile(dirName,'*.1.txt'));

filenames=struct2cell(pullfiles);

filematrix=transpose(filenames(1,:));
phvals(1,:)=cell2mat(filenames(1,:));

n=length(filematrix);
phvalscell={1,n};
for i=1:n-1;
    phvalscell{1,1}=phvals(1,1:3);
    phvalscell{1,i+1}=phvals(1,(i*9+1):(i*9+1)+2);
end

phvals1mat=(str2double(phvalscell))*0.01;

%Read in the uv data from each file into a matrix and then the wavelengths
%into another matrix

uvdata=zeros(1201,n);
for i=1:n;
    uvdata(:,i)=dlmread(fullfile(filepath,filematrix{i,1}),'\t',[20,2,1220,2]);
end

uvwavelength=dlmread(fullfile(filepath,filematrix{1,1}),'\t',[20,0,1220,0]);
%sizes of data matrices

[~,r]=size(uvdata);

%Subtract the first column from the last column (low pH from high pH)

uvdif=uvdata(:,r)-uvdata(:,1);
uvdifnorm=(uvdif)./(max(uvdif(1:1000)));

%Concatenate the UV wavelength with uvdif

uvplot=[uvwavelength,uvdifnorm];

[~,maxwaverow]=max(uvdif(1:1000,:));
maxwave=uvwavelength(maxwaverow,1);

[~,minwaverow]=min(uvdif(1:1000,:));
minwave=uvwavelength(minwaverow,1);

%plot the uv difference spectra
horzgridlinesx=200:320;
horzgridlinesy=zeros(121);

figure(6)
plot(horzgridlinesx,horzgridlinesy,':k','LineWidth',2)
hold on
legend off
hold on
plot(uvplot(:,1),uvplot(:,2),'k','LineWidth',2)
axis([220 320 -1 1.1])
set(get(gca,'XLabel'),'String','\lambda (nm)','FontName',...
    'Times New Roman','FontSize',32,'FontWeight','bold')
set(get(gca,'YLabel'),'String',...
    'A_{262} (a.u.)','FontName',...
    'Times New Roman','FontSize',32,'Fontweight','bold')
set(gca,'XTick',220:20:320)
set(gca,'XTickLabel',{'220','240','260','280','300','320'},...
    'FontSize',24,'FontName','Times New Roman')
set(gca,'YTick',-0.5:0.5:1)
set(gca,'YTickLabel',char('-0.5',' 0.0',' 0.5',' 1.0'),...
    'FontSize',24,'FontName','Times New Roman')

graphname9=strcat(filepath,strcat('/',sampleID,'_',...
    num2str(modelfuntoggle),'_UVdiffspectraLowRes'));
print(fullfile(graphname9),'-dtiff','-r200');

graphname10=strcat(filepath,strcat('/',sampleID,'_',...
    num2str(modelfuntoggle),'_UVdiffspectraHighRes'));
print(fullfile(graphname10),'-dtiff','-r1000');

fprintf('\n UV DIFF SPECTRA COMPLETE \n');

%print a results file with wavelength maximum and minimum as well as the pH
%values from which the difference spectra was compiled

resultspath=strcat(filepath,strcat('/',sampleID,'_',...
    num2str(modelfuntoggle),'_UVdiffResults.txt'));
fid=fopen(resultspath,'wt');

printsampleID=sampleID;
fprintf(fid,printsampleID);
fprintf(fid,'\n\n');
printmaxwave='max wavelength is %3.1f\n\n';
fprintf(fid,printmaxwave,maxwave);
printminwave='min wavelength is %3.1f\n\n';
fprintf(fid,printminwave,minwave);
pHmax='high pH spectra recorded at pH = %3.2f\n\n';
fprintf(fid,pHmax,phvals1mat(1,n));
pHmax='low pH spectra recorded at pH = %3.2f\n\n';
fprintf(fid,pHmax,phvals1mat(1,1));

fclose(fid);

fprintf('\n UVDIFF RESULTS COMPLETE \n');

%%%%%Print input file for timecourse
resultspath=strcat(filepath,strcat('/',sampleID,'_TimeCourseInput.txt'));
fid=fopen(resultspath,'wt');

printsampleID=sampleID;
fprintf(fid,printsampleID);
printcd='\nCD Analysis Wavelength\n%3.1f nm';
fprintf(fid,printcd,UVwavelength(row,1));
printrow='\nrownum\n%3.0f';
fprintf(fid,printrow,row);
printmax1='\n max1\n%3.15f';
printmin1='\n min1\n%3.15f';
printmax2='\n max2\n%3.15f';
printmin2='\n min2\n%3.15f';
printmax3='\n max3\n%3.15f';
printmin3='\n min3\n%3.15f';
fprintf(fid,printmax1,max1);
fprintf(fid,printmin1,min1);
fprintf(fid,printmax2,max2);
fprintf(fid,printmin2,min2);
fprintf(fid,printmax3,max3);
fprintf(fid,printmin3,min3);
fprintf(fid,'\nconcentration');
fprintf(fid,'\n%2.4f\t%1.4f',concdata);

fclose all;

%Updated Transition plot with biphasic regression
Ypred2calc=vertcat(0,Ypred2);
n1=length(Ypred2calc);
totalslope=sum(Ypred2(:,1));
sig=zeros(n1-1,1);
for i=3:n1
    sig(i-1,1)=sig(i-2,1)+Ypred2calc(i,1);
end
signorm=zeros(n1-1,1);
for i=1:n1-1
    signorm(i,1)=sig(i,1)/totalslope;
end
signormshiftandflip=-signorm+1;

Ypredupper2calc=signormshiftandflip+delta2;
Ypredlower2calc=signormshiftandflip-delta2;

figure(7)
box on
hold on
axis([5 8 -0.1 1.1])
%Plot the upper 95% CI with red dash and incorporate legend to denote it
plot(Testdata2,Ypredupper2calc,'--r','LineWidth',2,'Color','r')
hold on
legend({'\Delta_{95%}'},'FontName','Times New Roman','FontSize',24)
hold on
%Plot the regression line in black and the lower 95% CI
plot(Testdata2,signormshiftandflip,'-k','LineWidth',2,'Color','k')
hold on
plot(Testdata2,Ypredlower2calc,'--r','LineWidth',2,'Color','r')
hold on
%Plot the scatter data on top of the regressions
scatter(phtransitionplotALL(1:phx1,1),phtransitionplotALL(1:phx1,2),...
    35,'or','filled','MarkerEdgeColor','k','Linewidth',2)
hold on
scatter(phtransitionplotALL(phx1+1:(phx1+phx2),1),...
    phtransitionplotALL((phx1+1):(phx1+phx2),2),35,'or','filled',...
    'MarkerEdgeColor','k','Linewidth',2)
hold on
scatter(phtransitionplotALL((phx1+phx2+1):(phx1+phx2+phx3),1),...
    phtransitionplotALL((phx1+phx2+1):(phx1+phx2+phx3),2),35,'or',...
    'filled','MarkerEdgeColor','k','Linewidth',2)
hold on
%Set the axis labels and ticks
set(get(gca,'XLabel'),'String','pH','FontName',...
    'Times New Roman','FontSize',32,'FontWeight','bold')
set(get(gca,'YLabel'),'String',...
    '{\itx}_{iM}','FontName',...
    'Times New Roman','FontSize',32,'Fontweight','bold')
set(gca,'XTick',5:.5:8)
set(gca,'XTickLabel',char('5.0','5.5','6.0','6.5',...
    '7.0','7.5','8.0'),'FontSize',24,'FontName','Times New Roman')
set(gca,'YTick',0:0.2:1)
set(gca,'YTickLabel',char('0.0','0.2','0.4','0.6','0.8','1.0'),...
    'FontSize',24,'FontName','Times New Roman')
%Print the graph
graphname11=strcat(filepath,strcat('/',sampleID,'_',...
    num2str(modelfuntoggle),'_NumericalTransPlotLowRes'));
print(fullfile(graphname11),'-dtiff','-r200');

graphname12=strcat(filepath,strcat('/',sampleID,'_',...
    num2str(modelfuntoggle),'_NumericalTransPlotHighRes'));
print(fullfile(graphname12),'-dtiff','-r1000');

fprintf('\n NUMERICAL CONVERSION TO SIGMOID COMPLETE \n');

%Print a file with the Scatter Data, Regression Lines and 95% CI for use in
%a comparative plot for up and down data or more
resultspath=strcat(filepath,'/',sampleID,'_',num2str(modelfuntoggle),...
    '_GaussNumericalScatter1.txt');
fid=fopen(resultspath,'wt');


fprintf(fid,'%3.4f\t%3.4f\n',...
    [Testdata2,signormshiftandflip].');
fprintf(fid,'Upper 95%%CI Number of Pts and Data Set');
fprintf(fid,'\n81');
fprintf(fid,'\n%3.4f\t%3.4f',...
    [Testdata2,Ypredupper2calc].');
fprintf(fid,'\nLower 95%% CI Number of Pts and Data Set');
fprintf(fid,'81');
fprintf(fid,'\n%3.4f\t%3.4f',...
    [Testdata2,Ypredlower2calc].');

fclose all;
close all

fprintf('\n\n\n!...ANALYSIS COMPLETE...!\n\n\n\n');
end

