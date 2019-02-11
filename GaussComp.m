function [] = GaussComp(numfiles,parentdirectory,sampleID,statinterval)
%Plot comparisons of pHT graphs

%Determine the name of each file in the parent directory, the contents of
%which contain the numerical scatter output from PromptMasterCD

if numfiles==2
    f1=input('type the name of file 1\n');
    f2=input('type the name of file 2\n');
else if numfiles==3
        f1=input('type the name of file 1\n');
        f2=input('type the name of file 2\n');
        f3=input('type the name of file 3\n');
    else if numfiles==4
            f1=input('type the name of file 1\n');
            f2=input('type the name of file 2\n');
            f3=input('type the name of file 3\n');
            f4=input('type the name of file 4\n');
        else if numfiles==5
                f1=input('type the name of file 1\n');
                f2=input('type the name of file 2\n');
                f3=input('type the name of file 3\n');
                f4=input('type the name of file 4\n');
                f5=input('type the name of file 5\n');
            else if numfiles==6
                    f1=input('type the name of file 1\n');
                    f2=input('type the name of file 2\n');
                    f3=input('type the name of file 3\n');
                    f4=input('type the name of file 4\n');
                    f5=input('type the name of file 5\n');
                    f6=input('type the name of file (dashed) 6\n');
                else if numfiles==7
                        f1=input('type the name of file 1\n');
                        f2=input('type the name of file 2\n');
                        f3=input('type the name of file 3\n');
                        f4=input('type the name of file 4\n');
                        f5=input('type the name of file 5\n');
                        f6=input('type the name of file (dashed) 6\n');
                        f7=input('type the name of file (dashed) 7\n');
                    else if numfiles==8
                            f1=input('type the name of file 1\n');
                            f2=input('type the name of file 2\n');
                            f3=input('type the name of file 3\n');
                            f4=input('type the name of file 4\n');
                            f5=input('type the name of file 5\n');
                            f6=input('type the name of file 6\n');
                            f7=input('type the name of file 7\n');
                            f8=input('type the name of file 8\n');
                        end
                    end
                end
            end
        end
    end
end

%Construct the full filepaths for each numerical scatter output file

if numfiles==2
    f1path=struct2cell(dir(fullfile(strcat(parentdirectory,f1),...
        '*GaussPeakScatter.txt')));
    f2path=struct2cell(dir(fullfile(strcat(parentdirectory,f2),...
        '*GaussPeakScatter.txt')));
else if numfiles==3
        f1path=struct2cell(dir(fullfile(strcat(parentdirectory,f1),...
            '*GaussPeakScatter.txt')));
        f2path=struct2cell(dir(fullfile(strcat(parentdirectory,f2),...
            '*GaussPeakScatter.txt')));
        f3path=struct2cell(dir(fullfile(strcat(parentdirectory,f3),...
            '*GaussPeakScatter.txt')));
    else if numfiles==4
            f1path=struct2cell(dir(fullfile(strcat(parentdirectory,f1),...
                '*GaussPeakScatter.txt')));
            f2path=struct2cell(dir(fullfile(strcat(parentdirectory,f2),...
                '*GaussPeakScatter.txt')));
            f3path=struct2cell(dir(fullfile(strcat(parentdirectory,f3),...
                '*GaussPeakScatter.txt')));
            f4path=struct2cell(dir(fullfile(strcat(parentdirectory,f4),...
                '*GaussPeakScatter.txt')));
        else if numfiles==5
                f1path=struct2cell(dir(fullfile(...
                    strcat(parentdirectory,f1),'*GaussPeakScatter.txt')));
                f2path=struct2cell(dir(fullfile(...
                    strcat(parentdirectory,f2),'*GaussPeakScatter.txt')));
                f3path=struct2cell(dir(fullfile(...
                    strcat(parentdirectory,f3),'*GaussPeakScatter.txt')));
                f4path=struct2cell(dir(fullfile(...
                    strcat(parentdirectory,f4),'*GaussPeakScatter.txt')));
                f5path=struct2cell(dir(fullfile(...
                    strcat(parentdirectory,f5),'*GaussPeakScatter.txt')));
            else if numfiles==6
                    f1path=struct2cell(dir(fullfile(...
                        strcat(parentdirectory,f1),...
                        '*GaussPeakScatter.txt')));
                    f2path=struct2cell(dir(fullfile(...
                        strcat(parentdirectory,f2),...
                        '*GaussPeakScatter.txt')));
                    f3path=struct2cell(dir(fullfile(...
                        strcat(parentdirectory,f3),...
                        '*GaussPeakScatter.txt')));
                    f4path=struct2cell(dir(fullfile(...
                        strcat(parentdirectory,f4),...
                        '*GaussPeakScatter.txt')));
                    f5path=struct2cell(dir(fullfile(...
                        strcat(parentdirectory,f5),...
                        '*GaussPeakScatter.txt')));
                    f6path=struct2cell(dir(fullfile(...
                        strcat(parentdirectory,f6),...
                        '*GaussPeakScatter.txt')));
                else if numfiles==7
                        f1path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f1),...
                            '*GaussPeakScatter.txt')));
                        f2path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f2),...
                            '*GaussPeakScatter.txt')));
                        f3path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f3),...
                            '*GaussPeakScatter.txt')));
                        f4path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f4),...
                            '*GaussPeakScatter.txt')));
                        f5path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f5),...
                            '*GaussPeakScatter.txt')));
                        f6path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f6),...
                            '*GaussPeakScatter.txt')));
                        f7path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f7),...
                            '*GaussPeakScatter.txt')));
                    else if numfiles==8
                            f1path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f1),...
                            '*GaussPeakScatter.txt')));
                        f2path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f2),...
                            '*GaussPeakScatter.txt')));
                        f3path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f3),...
                            '*GaussPeakScatter.txt')));
                        f4path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f4),...
                            '*GaussPeakScatter.txt')));
                        f5path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f5),...
                            '*GaussPeakScatter.txt')));
                        f6path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f6),...
                            '*GaussPeakScatter.txt')));
                        f7path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f7),...
                            '*GaussPeakScatter.txt')));
                        f8path=struct2cell(dir(fullfile(...
                            strcat(parentdirectory,f8),...
                            '*GaussPeakScatter.txt')));
                        end
                    end
                end
            end
        end
    end
end

%Read in the data for the Comparison plot without the data points

if numfiles==2
    f1data=dlmread(fullfile(parentdirectory,...
        f1,f1path{1,1}),'\t',[1,0,81,3]);
    f2data=dlmread(fullfile(parentdirectory,...
        f2,f2path{1,1}),'\t',[1,0,81,3]);
else if numfiles==3
        f1data=dlmread(fullfile(parentdirectory,...
            f1,f1path{1,1}),'\t',[1,0,81,3]);
        f2data=dlmread(fullfile(parentdirectory,...
            f2,f2path{1,1}),'\t',[1,0,81,3]);
        f3data=dlmread(fullfile(parentdirectory,...
            f3,f3path{1,1}),'\t',[1,0,81,3]);
    else if numfiles==4
            f1data=dlmread(fullfile(parentdirectory,...
                f1,f1path{1,1}),'\t',[1,0,81,3]);
            f2data=dlmread(fullfile(parentdirectory,...
                f2,f2path{1,1}),'\t',[1,0,81,3]);
            f3data=dlmread(fullfile(parentdirectory,...
                f3,f3path{1,1}),'\t',[1,0,81,3]);
            f4data=dlmread(fullfile(parentdirectory,...
                f4,f4path{1,1}),'\t',[1,0,81,3]);
        else if numfiles==5
                f1data=dlmread(fullfile(parentdirectory,...
                    f1,f1path{1,1}),'\t',[1,0,81,3]);
                f2data=dlmread(fullfile(parentdirectory,...
                    f2,f2path{1,1}),'\t',[1,0,81,3]);
                f3data=dlmread(fullfile(parentdirectory,...
                    f3,f3path{1,1}),'\t',[1,0,81,3]);
                f4data=dlmread(fullfile(parentdirectory,...
                    f4,f4path{1,1}),'\t',[1,0,81,3]);
                f5data=dlmread(fullfile(parentdirectory,...
                    f5,f5path{1,1}),'\t',[1,0,81,3]);
            else if numfiles==6
                    f1data=dlmread(fullfile(parentdirectory,...
                        f1,f1path{1,1}),'\t',[1,0,81,3]);
                    f2data=dlmread(fullfile(parentdirectory,...
                        f2,f2path{1,1}),'\t',[1,0,81,3]);
                    f3data=dlmread(fullfile(parentdirectory,...
                        f3,f3path{1,1}),'\t',[1,0,81,3]);
                    f4data=dlmread(fullfile(parentdirectory,...
                        f4,f4path{1,1}),'\t',[1,0,81,3]);
                    f5data=dlmread(fullfile(parentdirectory,...
                        f5,f5path{1,1}),'\t',[1,0,81,3]);
                    f6data=dlmread(fullfile(parentdirectory,...
                        f6,f6path{1,1}),'\t',[1,0,81,3]);
                else if numfiles==7
                        f1data=dlmread(fullfile(parentdirectory,...
                            f1,f1path{1,1}),'\t',[1,0,81,3]);
                        f2data=dlmread(fullfile(parentdirectory,...
                            f2,f2path{1,1}),'\t',[1,0,81,3]);
                        f3data=dlmread(fullfile(parentdirectory,...
                            f3,f3path{1,1}),'\t',[1,0,81,3]);
                        f4data=dlmread(fullfile(parentdirectory,...
                            f4,f4path{1,1}),'\t',[1,0,81,3]);
                        f5data=dlmread(fullfile(parentdirectory,...
                            f5,f5path{1,1}),'\t',[1,0,81,3]);
                        f6data=dlmread(fullfile(parentdirectory,...
                            f6,f6path{1,1}),'\t',[1,0,81,3]);
                        f7data=dlmread(fullfile(parentdirectory,...
                            f7,f7path{1,1}),'\t',[1,0,81,3]);
                    else if numfiles==8
                            f1data=dlmread(fullfile(parentdirectory,...
                            f1,f1path{1,1}),'\t',[1,0,81,3]);
                        f2data=dlmread(fullfile(parentdirectory,...
                            f2,f2path{1,1}),'\t',[1,0,81,3]);
                        f3data=dlmread(fullfile(parentdirectory,...
                            f3,f3path{1,1}),'\t',[1,0,81,3]);
                        f4data=dlmread(fullfile(parentdirectory,...
                            f4,f4path{1,1}),'\t',[1,0,81,3]);
                        f5data=dlmread(fullfile(parentdirectory,...
                            f5,f5path{1,1}),'\t',[1,0,81,3]);
                        f6data=dlmread(fullfile(parentdirectory,...
                            f6,f6path{1,1}),'\t',[1,0,81,3]);
                        f7data=dlmread(fullfile(parentdirectory,...
                            f7,f7path{1,1}),'\t',[1,0,81,3]);
                        f8data=dlmread(fullfile(parentdirectory,...
                            f8,f8path{1,1}),'\t',[1,0,81,3]);
                        end
                    end
                end
            end
        end
    end
end

%Plot the data for visual comparison of the regression intervals

figure(1)
box on
hold on
axis([5 8 -0.1 1.1])
hold on

if numfiles==2
    plot(f1data(:,1),f1data(:,2),'-ok','LineWidth',3)
    hold on
    plot(f2data(:,1),f2data(:,2),'-or','LineWidth',3)
    hold on
else if numfiles==3
        plot(f1data(:,1),f1data(:,2),'-k','LineWidth',3)
        hold on
        plot(f2data(:,1),f2data(:,2),'-r','LineWidth',3)
        hold on
        plot(f3data(:,1),f3data(:,2),'-b','LineWidth',3)
        hold on
        if statinterval==1
            plot(f1data(:,1),f1data(:,3),'--k','LineWidth',2)
            plot(f1data(:,1),f1data(:,4),'--k','LineWidth',2)
            plot(f2data(:,1),f2data(:,3),'--r','LineWidth',2)
            plot(f2data(:,1),f2data(:,4),'--r','LineWidth',2)
            plot(f3data(:,1),f3data(:,3),'--b','LineWidth',2)
            plot(f3data(:,1),f3data(:,4),'--b','LineWidth',2)
        else if statinterval==0
            end
        end
    else if numfiles==4
            plot(f1data(:,1),f1data(:,2),'-k','LineWidth',2)
            hold on
            plot(f2data(:,1),f2data(:,2),'-b','LineWidth',2)
            hold on
            plot(f3data(:,1),f3data(:,2),'-','Color',[0 0.5 0],'LineWidth',2)
            hold on
            plot(f4data(:,1),f4data(:,2),'-r','LineWidth',2)
            hold on
        else if numfiles==5
                plot(f1data(:,1),f1data(:,2),'k-','LineWidth',2)
                hold on
                plot(f2data(:,1),f2data(:,2),'r-','LineWidth',2)
                hold on
                plot(f3data(:,1),f3data(:,2),'b-','LineWidth',2)
                hold on
                plot(f4data(:,1),f4data(:,2),'-','Color',[0.91 0.41 0.17],'LineWidth',2)
                hold on
                plot(f5data(:,1),f5data(:,2),'m-','LineWidth',2)
                hold on
            else if numfiles==6
                    plot(f1data(:,1),f1data(:,2),'-or','LineWidth',2)
                    hold on
                    plot(f2data(:,1),f2data(:,2),'-og','LineWidth',2)
                    hold on
                    plot(f3data(:,1),f3data(:,2),'-ob','LineWidth',2)
                    hold on
                    plot(f4data(:,1),f4data(:,2),'-xb','LineWidth',2)
                    hold on
                    plot(f5data(:,1),f5data(:,2),'-xg','LineWidth',2)
                    hold on
                    plot(f6data(:,1),f6data(:,2),'-ok','LineWidth',2)
                    hold on
                else if numfiles==7
                        plot(f1data(:,1),f1data(:,2),'k-','LineWidth',2)
                        hold on
                        plot(f2data(:,1),f2data(:,2),'r-','LineWidth',2)
                        hold on
                        plot(f3data(:,1),f3data(:,2),'b-','LineWidth',2)
                        hold on
                        plot(f4data(:,1),f4data(:,2),'-','Color',[0.91 0.41 0.17],'LineWidth',2)
                        hold on
                        plot(f5data(:,1),f5data(:,2),'m-','LineWidth',2)
                        hold on
                        plot(f6data(:,1),f6data(:,2),'k-.','LineWidth',1)
                        hold on
                        plot(f7data(:,1),f7data(:,2),'k--','LineWidth',1)
                        hold on
                    else if numfiles==8
                            plot(f1data(:,1),f1data(:,2),'b-','LineWidth',2)
                            hold on
                            plot(f2data(:,1),f2data(:,2),'-','Color',[0 0.5 0],'LineWidth',2)
                            hold on
                            plot(f3data(:,1),f3data(:,2),'-..','Color',[0 0.5 0],'LineWidth',2)
                            hold on
                            plot(f4data(:,1),f4data(:,2),'m-','LineWidth',2)
                            hold on
                            plot(f5data(:,1),f5data(:,2),'m-..','LineWidth',2)
                            hold on
                            plot(f6data(:,1),f6data(:,2),'-','Color',[0.5 0 0.8],'LineWidth',2)
                            hold on
                            plot(f7data(:,1),f7data(:,2),'k-','LineWidth',3)
                            hold on
                            plot(f8data(:,1),f8data(:,2),'r-','LineWidth',3)
                            hold on
                        end
                    end
                end
            end
        end
    end
end

hold on
set(get(gca,'XLabel'),'String','pH','FontName',...
    'Times New Roman','FontSize',48,'FontWeight','bold')
set(get(gca,'YLabel'),'String','\Delta {\it{x}}/\Delta pH',...
    'FontName','Times New Roman','FontSize',48,'FontWeight','bold')
set(gca,'XTick',5:1:8)
set(gca,'XTickLabel',char('5.0','6.0','7.0','8.0'),...
    'FontSize',36,'FontName','Times New Roman')
set(gca,'YTick',0:0.25:1)
set(gca,'YTickLabel',char('0.00','0.25','0.50','0.75','1.00'),...
    'FontSize',36,'FontName','Times New Roman')

graphname1=strcat(parentdirectory,strcat(sampleID,'_GPLowRes'));
print(fullfile(graphname1),'-dtiff','-r200');

graphname1=strcat(parentdirectory,strcat(sampleID,'_GPHighRes'));
print(fullfile(graphname1),'-dtiff','-r1000');

end

