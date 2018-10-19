%script to read example data product

clear all; close all;

addpath E:\readSWOTBeta     %replace this line with local path to : https://github.com/mikedurand/readSWOTBeta
addpath E:\SWOTAprimeCalcs  %replace this line with local path to : https://github.com/mikedurand/SWOTAprimeCalcs

DataDir='C:\Users\wei.263\Desktop\L2_HR_Final\L2_HR_Sac_data_product_6 - Copy\';
%DataDir='E:\SWOTSAC1051_NEWellip\RO\Sac_v7\';
%DataDir='E:\SWOTPO1051\Po_floodplain_study\RiverObs\Po_v2\'L2_HR_Final\L2_HR_Sac_data_product_6;  %replace with path to example dataset: go.osu.edu/swotbeta
River='Sac';

% define passes and cycle variables
DefinePassCycle_2;

% read reaches & nodes
[SWOTReaches,TrueReaches,SWOTNodes,TrueNodes]=ReadShapeData(Cycles,Passes,Dates,DataDir,River);

% extract reaches & nodes
ExtractHWS;

% error stats
CalcErrStats_2;

% check partial observation
for i=1:nObs
    for j=1:nReach
        fracObs(j,i)=sum(~isnan(Hn(NodesInReach{j},i)))/sum(NodesInReach{j});
    end
end

% for i=1:nPass
%     for j=1:nReach
%         fracObsPass(j,i)=median(fracObs(j,[i:3:nObs]),2);
%     end
% end

for i=1:nPass
    for j=1:nReach
        fracObsPass(j,i)=median(fracObs(j,ind.(['Pass' num2str(Pass(i))])),2);
    end
end

zone = utmzone(latn(1),lonn(1));
utmstruct=defaultm('utm');
utmstruct.zone=zone;
utmstruct=defaultm(utmstruct);
[X, Y] = mfwdtran(utmstruct,latn(:,1), lonn(:,1));
FD = GetFlowDist(nan(size(X)),Y,-9999);
for i=1:nObs
    [~,iSort{i}]=sort(ReachIDnt(:,i));
end

% for i=1:nReach
%     FDr(i,1)=FD(find(ReachIDnt(:,1)==i,1,'first'));
%     FDr(i,2)=FD(find(ReachIDnt(:,1)==i,1,'last'));    
%     RL(i)=FDr(i,2)-FDr(i,1);
% end

for i=1:length(ObsReachIDs)
    FDr(i,1)=FD(find(ReachIDnt_2(:,1)==ObsReachIDs(i),1,'first'));
    FDr(i,2)=FD(find(ReachIDnt_2(:,1)==ObsReachIDs(i),1,'last'));  
    RL(i)=FDr(i,2)-FDr(i,1);
end

if ObsReachIDs(1) ~= 1
    RL = RL(ObsReachIDs(1):ObsReachIDs(end));
end
%% overview plots
figure(1)
for i =1:nPass
   firstDay(i) = find(ind.(['Pass' num2str(Pass(i))]), 1, 'first');
end
for i = 1:nPass
    subplot( 1 ,nPass, i)
    mapshow(TrueReaches(firstDay(i)).S,'Color','Blue')
    mapshow(SWOTNodes(firstDay(i)).S,'Color','Red','Marker','o')
    set(gca,'FontSize',14)
    title(['Pass' num2str(Pass(i))])
end

figure(2) %missing data and reach overveiw
subplot(311)
pcolor(Hr)
colorbar
set(gca,'FontSize',14)
title('Reach height [m]')
ylabel('Reach #')
subplot(312)
pcolor(Sr)
set(gca,'FontSize',14)
title('Reach slope [mm/km]')
ylabel('Reach #')
colorbar
subplot(313)
pcolor(Wr)
colorbar
set(gca,'FontSize',14)
title('Reach width [m]')
xlabel('Observation #')
ylabel('Reach #')

figure(3)
subplot(131)
boxplot(Err.RchHeight.All.Allv)
set(gca,'FontSize',14)
title('Height')
ylabel('Height error, m')
subplot(132)
boxplot(Err.RchWidth.All.Allv)
set(gca,'FontSize',14)
title('Width')
ylabel('Width error, m')
subplot(133)
boxplot(Err.RchSlope.All.Allv)
set(gca,'FontSize',14)
title('Slope')
ylabel('Slope error, mm/km')

figure(30)
subplot(121)
boxplot(Err.RchSlope.All.Allv)
title(['RMSE,' num2str(round(Err.RchSlope.All.RMSE,2))])
set(gca,'FontSize',14)
ylabel('Slope error, mm/km')
subplot(122)
boxplot(Err.RchSlopeEn.All.Allv)
title(['RMSE,' num2str(round(Err.RchSlopeEn.All.RMSE,2))])
set(gca,'FontSize',14)
ylabel('Enhanced Slope error, mm/km')

figure(4)
plot(1:nReach,fracObsPass,'o-','LineWidth',2)
set(gca,'FontSize',14)
xlabel('Reach #')
ylabel('Fraction Observed (Median across cycles)')
legend(['Pass' num2str(Pass(1))],['Pass' num2str(Pass(2))],['Pass' num2str(Pass(3))],'Location','Best')
grid on


%% height plots

figure(5) %reach height error by pass
for i = 1:nPass
    subplot(1, nPass, i)
    boxplot(Err.RchHeight.(['Pass' num2str(Pass(i))]).Allv)
    set(gca,'FontSize',14)
    ylabel('Height error, m')
    title(['Pass' num2str(Pass(i))])

end

figure(6) %reach height errors by cycle and node
for i = 1:nPass
    subplot(nPass, 1, i)
    if size(Err.RchHeight.(['Pass' num2str(Pass(i))]).All, 2) == 1 
        pcolor([Err.RchHeight.(['Pass' num2str(Pass(i))]).All, ...
                Err.RchHeight.(['Pass' num2str(Pass(i))]).All])
    else
        pcolor(Err.RchHeight.(['Pass' num2str(Pass(i))]).All)
    end
    colorbar
    set(gca,'CLim',[-0.2 0.2])
    xlabel('cycle')
    ylabel('reach')
    title(['Height errors Pass' num2str(Pass(i))])
end

figure(7) %it's not cross-track distance!
for i =1:nPass
    plot(XTDr(:,i)./1000,rms(Err.RchHeight.(['Pass' num2str(Pass(i))]).All,2),'o','LineWidth',2); 
    hold on;
end
hold off
set(gca,'FontSize',14)
legend(['Pass' num2str(Pass(1))],['Pass' num2str(Pass(2))],['Pass' num2str(Pass(3))],'Location','Best')
xlabel('Reach-averaged cross-track distance, km')
ylabel('Height RMS error')
grid on

%% partially-observed reach sample calculation
figure(8)
i8=find(NodesInReach{8});
i10=find(NodesInReach{10});
subplot(2,1,1)
plot(FD(i8),Hnt(i8,1),FD(i8),Hn(i8,1),'LineWidth',2)
set(gca,'FontSize',14)
ylabel('Water surface elevation, m')
title('Reach 8, Pass 553, Cycle 1')
legend('True','SWOT')
grid on;
subplot(2,1,2)
plot(FD(i10),Hnt(i10,1),FD(i10),Hn(i10,1),'LineWidth',2)
% plot(NodeID(i5),Hnt(i5,3),NodeID(i5),Hn(i5,3),'LineWidth',2)
set(gca,'FontSize',14)
xlabel('Flow distance, s, km')
% xlabel('Node ID')
ylabel('Water surface elevation, m')
title('Reach 10, Pass 553, Cycle 1')
legend('True','SWOT')
grid on;

% i1=find(NodesInReach{1});
% x_t = FD(i5)-51.8;
% x_t2 = x_t - max(x_t)/2;
% cef_t = polyfit(x_t2,Hnt(i5,9)',1)
% x=FD(i5)-51.8;
% x2=x-max(x)/2;
% y=Hn(i5,9)';
% idx=~isnan (Hn(i5,9));
% cef_est = polyfit(x2(idx),y(idx),1)
% 
% Method='Current';
% [TrueSlope,TrueHeight] = SampleCalcReachAvg(FD(i5).*1000,Hnt(i5,3)',Method);
% [SWOTSlope,SWOTHeight] = SampleCalcReachAvg(FD(i5).*1000,Hn(i5,3)',Method);
% 
% disp(['RiverObs vs Sample calc: truth: ' num2str(TrueReaches(3).A(5).height) ...
%     ' vs ' num2str(TrueHeight)])
% disp(['RiverObs vs Sample calc: SWOT: ' num2str(SWOTReaches(3).A(5).height) ...
%     ' vs ' num2str(SWOTHeight)])
% 
% Method='Fixed';
% [SWOTSlopeFix,SWOTHeightFix] = SampleCalcReachAvg(FD(i5).*1000,Hn(i5,3)',Method);
% 
% disp(['RiverObs True vs SWOT fixed Sample calc: ' num2str(TrueReaches(3).A(5).height) ...
%     ' vs ' num2str(SWOTHeightFix)])


figure(9) %reach slope errors by cycle and node
for i =1:nPass
    subplot(nPass, 1, i)
    if size(Err.RchSlope.(['Pass' num2str(Pass(i))]).All, 2) == 1 
        pcolor([Err.RchSlope.(['Pass' num2str(Pass(i))]).All, ...
                Err.RchSlope.(['Pass' num2str(Pass(i))]).All])
    else
        pcolor(Err.RchSlope.(['Pass' num2str(Pass(i))]).All)
    end
    colorbar
    set(gca,'CLim',[-10 10])
    xlabel('cycle')
    ylabel('reach')
    title(['Slope errors: Pass' num2str(Pass(i))])
end

figure(10)
for i =1:nPass
    subplot(1, nPass, i)
    boxplot(Err.RchSlope.(['Pass' num2str(Pass(i))]).Allv)
    set(gca,'FontSize',12)
    title(['RMSE,' num2str(round(Err.RchSlope.(['Pass' num2str(Pass(i))]).RMSE,2))])
    ylabel('Slope error, mm/km')
end

figure(11)
for i =1:nPass
    subplot(1, nPass, i)
    boxplot(Err.RchSlopeEn.(['Pass' num2str(Pass(i))]).Allv)
    set(gca,'FontSize',12)
    ylabel('Enhanced Slope error, mm/km')
    title(['RMSE,' num2str(round(Err.RchSlopeEn.(['Pass' num2str(Pass(i))]).RMSE,2))])
end

figure(12) %reach slope errors by cycle and node
for i =1:nPass
    subplot(nPass, 1, i)
    
    if size(Err.RchSlope.(['Pass' num2str(Pass(i))]).All, 2) == 1 
        pcolor([Err.RchWidth.(['Pass' num2str(Pass(i))]).All, ...
                Err.RchWidth.(['Pass' num2str(Pass(i))]).All])
    else
        pcolor(Err.RchWidth.(['Pass' num2str(Pass(i))]).All)
    end
    colorbar
    % set(gca,'CLim',[-10 10])
    xlabel('cycle')
    ylabel('reach')
    title(['Width errors: Pass' num2str(Pass(i))])
end

figure(13)
[f13_x,ind_x]=sort(cell2mat(Dates));
plot(dates(ind_x),Hrt([6 7],ind_x),'b-',dates(ind_x),Hr([6 7],ind_x),'ro','MarkerSize',8)
set(gca,'FontSize',14)
datetick
ylabel('Water surface elevation, m')
grid

figure(14)
subplot(211)
plot(dates(ind_x),Srt([6 7],ind_x),'b-',dates(ind_x),Sr([6 7],ind_x),'ro','MarkerSize',8)
set(gca,'FontSize',14)
datetick
ylabel('Slope, mm/km')
title('Reach 6 and 7')
grid
subplot(212)
plot(dates(ind_x),Srte([6 7],ind_x),'b-',dates(ind_x),Sre([6 7],ind_x),'ro','MarkerSize',8)
set(gca,'FontSize',14)
datetick
ylabel('Enhanced slope')
grid

figure(15)
plot(dates(ind_x),Wrt([6,7],ind_x),'b-',dates(ind_x),Wr([6,7],ind_x),'ro','MarkerSize',8)
set(gca,'FontSize',14)
datetick
ylabel('Width, m')
title('Reach 6 and 7')
grid

figure(16)
subplot(221)
r=3;
plot(Hr(r,:),Wr(r,:),'o','LineWidth',2)
set(gca,'FontSize',14)
xlabel('Height, m'); ylabel('Width, m')
title('Reach 3')
grid on;
subplot(222)
r=4;
plot(Hr(r,:),Wr(r,:),'o','LineWidth',2)
set(gca,'FontSize',14)
xlabel('Height, m'); ylabel('Width, m')
title('Reach 4')
grid on;
subplot(223)
r=5;
plot(Hr(r,:),Wr(r,:),'o','LineWidth',2)
set(gca,'FontSize',14)
xlabel('Height, m'); ylabel('Width, m')
title('Reach 5')
grid on;
subplot(224)
r=10;
plot(Hr(r,:),Wr(r,:),'o','LineWidth',2)
set(gca,'FontSize',14)
xlabel('Height, m'); ylabel('Width, m')
title('Reach 10')
grid on;

for i=1:nReach
    for j=1:nPass
        RMSESlope(i).(['Pass' num2str(Pass(j))])=(nanmean ...
                      (Err.RchSlope.(['Pass' num2str(Pass(j))]).All(i,:)).^2+ ...
                       nanstd(Err.RchSlope.(['Pass' num2str(Pass(j))]).All(i,:)).^2).^.5;
                   
        RMSEHeight(i).(['Pass' num2str(Pass(j))])=(nanmean ...
                      (Err.RchSlope.(['Pass' num2str(Pass(j))]).All(i,:)).^2+ ...
                       nanstd(Err.RchHeight.(['Pass' num2str(Pass(j))]).All(i,:)).^2).^.5; 
    end
end

figure(31)
plot(RL,[RMSESlope.(['Pass' num2str(Pass(1))])],'o', ...
     RL,[RMSESlope.(['Pass' num2str(Pass(2))])],'o',...
     RL,[RMSESlope.(['Pass' num2str(Pass(3))])],'o','LineWidth',2) 
 
%% nodes
figure(17)
subplot(121)
cdfplot(nPix(:))
subplot(122)
plot(nPix(:),Wn(:)-Wnt(:),'.')
set(gca,'FontSize',14)
grid on
xlabel('Number of pixels used to compute inundated area')
ylabel('Width error, m')

figure(18) 
plot(FD,Hn(:,1),'.'); hold on;
plot(FD,Hn(:,2),'.');
plot(FD,Hn(:,3),'.'); hold off;
set(gca,'FontSize',14)
xlabel('Flow distance, km')
ylabel('Water surface elevation, m')
legend(datestr(dates(1:3)'),'Location','Best')
title(['Passes ' num2str(Pass(1)) ' , ' ...
       num2str(Pass(2)) ' and ' num2str(Pass(3)), ' from cycle 1'])

figure(19)
hist(Err.NodeHeight.Allv,50)
set(gca,'FontSize',14)
ylabel('Height Error, m')
title('SWOT height errors for nodes')

figure(20)
plot(FD,Err.NodeHeight.All,'.')
set(gca,'FontSize',14)
xlabel('Flow distance, km')
ylabel('Height Error, m')
title('SWOT height errors for nodes')

figure(21)
for i=1:nPass
    subplot(1, nPass, i)
    xt=XTDn(:,ind.(['Pass' num2str(Pass(i))]));
    hist(xt(:)./1000)
    set(gca,'FontSize',14)
    title(['Pass' num2str(Pass(i))])
    ylabel('Count')
    xlabel('Cross-track dist., km')
end

figure(22)
C=get(groot,'defaultAxesColorOrder');
plot(FD,Wnt(:,[1 2 3]),'LineWidth',2); hold on;
h=plot(FD,Wn(:,[1 2 3]),'.','LineWidth',2); hold off;
for i=1:3
    set(h(i),'Color',C(i,:))
end
set(gca,'FontSize',14)
xlabel('Flow distance, km')
ylabel('Width, m')
title('All Passes from Cycle 1')
legend(datestr(dates(1:3)'),'Location','Best')

figure(23)
XTDnBinBound=[10:62].*1000;
for i=1:length(XTDnBinBound)-1
    j=XTDn(:)>=XTDnBinBound(i) & XTDn(:)<XTDnBinBound(i+1) & ~isnan(Err.NodeWidth.All(:));
    RmsNodeWidthXT(i)=rms(Err.NodeWidth.All(j)); 
end
plot(XTDnBinBound(1:end-1)./1000+.5,RmsNodeWidthXT,'LineWidth',2)
set(gca,'FontSize',14)
xlabel('Cross-track distance, km')
ylabel('Node Width RMSE, m')
title('Node widths')
grid on;

figure(24)
XTDnBinBound=[10:62].*1000;
for i=1:length(XTDnBinBound)-1
    j=XTDn(:)>=XTDnBinBound(i) & XTDn(:)<XTDnBinBound(i+1) & ~isnan(Err.NodeHeight.All(:));
    RmsNodeHeightXT(i)=rms(Err.NodeHeight.All(j)); 
end
plot(XTDnBinBound(1:end-1)./1000+.5,RmsNodeHeightXT.*100,'o-','LineWidth',2)
set(gca,'FontSize',14)
xlabel('Cross-track distance, km')
ylabel('Node height RMSE, cm')
grid on

figure(25)
i=3;  
% j=2:nNode-1;
j=1:nNode;
plot(FD(j),Wnt(j,i),FD(j),Wn(j,i))
set(gca,'FontSize',14)
legend('True','SWOT','Location','Best')
xlabel('Flow distance, km')
ylabel('Node widths, m')
title(['Pass ' num2str(Passes(i)) ', Cycle ' num2str(Cycles(i))])

figure(26)
wcut=100;
hist(Err.NodeWidth.Allv(abs(Err.NodeWidth.Allv)<wcut),50)

figure(27)
hist(Wnt(:),50)
set(gca,'FontSize',14)
xlabel('River Width, m')
ylabel('Count')
title('True Node Widths')

figure(28)
nids=[100 101];
plot(dates,Wnt(nids,:),'x-','LineWidth',2); hold on;
h=plot(dates,Wn(nids,:),'o','LineWidth',2); hold off;
for i=1:2
    set(h(i),'Color',C(i,:))
end
set(gca,'FontSize',14)
datetick
ylabel('River width, m')
legend('True node #100','True node #101','SWOT node #100','SWOT node #101','Location','Best')

%% et al.: le geoid!
figure(29)
plot(FD,[TrueNodes(1).A.geoid_hght])