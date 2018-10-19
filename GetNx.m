%This script figures out all unique reach & node ids in the database, and
%indices the rows of the reach or node data structure, to allow for
%comparison across passes without all data present. i.e. it handles
%unobserved reaches and nodes

function [NxReach,ObsReachIDs,NxNode,ObsNodeIDs] = GetNx(Reaches,Nodes)

npass=length(Reaches);             %number of passes / obs times

ReachIDs=[];
for i=1:npass
    ReachIDs=[ReachIDs; [Reaches(i).A.reach_id]';];
end

ObsReachIDs=unique(ReachIDs);
NxReach=length(ObsReachIDs);

NodeIDs=[];
for i=1:npass
    NodeIDs=[NodeIDs; [[Nodes(i).A.reach_id]' [Nodes(i).A.node_id]'];];            
end

ObsNodeIDs=unique(NodeIDs,'rows');
NxNode=length(ObsNodeIDs);


return
