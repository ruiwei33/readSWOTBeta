
[NxReach,ObsReachIDs,NxNode,ObsNodeIDs] = GetNx(TrueReaches,TrueNodes);

Hrt=ExtractData(TrueReaches,'height',NODATA,NxReach,ObsReachIDs,'Reach');
Hr=ExtractData(SWOTReaches,'height',NODATA,NxReach,ObsReachIDs,'Reach');
Srt=ExtractData(TrueReaches,'slope2',NODATA,NxReach,ObsReachIDs,'Reach'); %slope2
Sr=ExtractData(SWOTReaches,'slope',NODATA,NxReach,ObsReachIDs,'Reach');
Sre=ExtractData(SWOTReaches,'slope2',NODATA,NxReach,ObsReachIDs,'Reach');
Srte=ExtractData(TrueReaches,'slope2',NODATA,NxReach,ObsReachIDs,'Reach');
Wrt=ExtractData(TrueReaches,'width',NODATA,NxReach,ObsReachIDs,'Reach');
Wr=ExtractData(SWOTReaches,'width',NODATA,NxReach,ObsReachIDs,'Reach');
% set reach width with <50% obs to nan, didn't to do this in RiverObs 
idxr_nan = isnan(Hr);
Wr(idxr_nan) = NaN;

%dAr=ExtractData(SWOTReaches,'del_X_Area',NODATA,NxReach,ObsReachIDs,'Reach');

Hnt=ExtractData(TrueNodes,'height',NODATA,NxNode,ObsNodeIDs,'Node');%+ExtractData(TrueNodes,'geoid_hght',NODATA,NxNode,ObsNodeIDs,'Node');
Hn=ExtractData(SWOTNodes,'height',NODATA,NxNode,ObsNodeIDs,'Node');%+ExtractData(SWOTNodes,'geoid_hght',NODATA,NxNode,ObsNodeIDs,'Node');
Wnt=ExtractData(TrueNodes,'width',NODATA,NxNode,ObsNodeIDs,'Node');
An=ExtractData(SWOTNodes,'area_detct',NODATA,NxNode,ObsNodeIDs,'Node');
Ant=ExtractData(TrueNodes,'area_detct',NODATA,NxNode,ObsNodeIDs,'Node');
Wn=ExtractData(SWOTNodes,'width',NODATA,NxNode,ObsNodeIDs,'Node');
XTDn=ExtractData(SWOTNodes,'xtrk_dist',NODATA,NxNode,ObsNodeIDs,'Node');
latn=ExtractData(TrueNodes,'latitude',NODATA,NxNode,ObsNodeIDs,'Node');
lonn=ExtractData(TrueNodes,'longitude',NODATA,NxNode,ObsNodeIDs,'Node');
ReachIDnt=ExtractData(TrueNodes,'reach_id',NODATA,NxNode,ObsNodeIDs,'Node');
NodeID=ExtractData(TrueNodes,'node_id',NODATA,NxNode,ObsNodeIDs,'Node');
nPix=ExtractData(SWOTNodes,'n_good_pix',NODATA,NxNode,ObsNodeIDs,'Node');

% X=ExtractData(SWOTNodes,'X',NODATA,NxNode,ObsNodeIDs,'Node');
% Y=ExtractData(SWOTNodes,'Y',NODATA,NxNode,ObsNodeIDs,'Node');
% Xt=ExtractData(TrueNodes,'X',NODATA,NxNode,ObsNodeIDs,'Node');
% Yt=ExtractData(TrueNodes,'Y',NODATA,NxNode,ObsNodeIDs,'Node');
% set node width with <100 obs to nan, didn't to do this in RiverObs 
idxn_nan = isnan(Hn);
Wn(idxn_nan) = NaN;

nReach=size(Hr,1);
nNode=size(Hn,1);

for i=1:size(ReachIDnt,1)
    row=[];
    row=ReachIDnt(i,:);
    nonnan=[];
    nonnan=row(find(~isnan(row)));
    ReachIDnt_2(i,1)=nonnan(1);
end

% for i=1:nReach % mike
%     NodesInReach{i}=[TrueNodes(1).A.reach_id]==i;
%     nNodesInReach(i)=sum(NodesInReach{i});
%     for j=1:nPass
%         XTDr(i,j)=nanmean(XTDn(NodesInReach{i},j));
%     end
% end

for i=1:nReach   % rui
    NodesInReach{i}=ReachIDnt_2==ObsReachIDs(i); 
    nNodesInReach(i)=sum(NodesInReach{i});
    for j=1:nPass
        XTDr(i,j)=nanmean(XTDn(NodesInReach{i},j));
    end
end

% RL=ExtractData(SWOTReaches,'Length',NODATA,NxReach,ObsReachIDs,'Reach'); 
% RL=RL(:,1);
% FD = GetFlowDist(xn,yn,NODATA);