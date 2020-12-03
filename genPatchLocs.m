% Generate the edge locations for a combination of `normal'
% patches, of width pWid, and `shock' patches of width spWid.
% The total number of patches is (at least) nPatch. The shock locations,
% i.e. midpoints for shock patches, are input in sLocs, and the
% domain boundaries are in xBnds.
% 
% Outputs:
% edges: location of patch edges in an nPatch x 2 array
% macroNodes: locations of the macroscale points. For normal
% patches, this is the centre of the patch. For shock patches,
% there are two nodes a distance spMacroDist from the edges.
% 
% sNodeInd: logical vector, same length as macroNodes. 0 if the node is
% associated with a normal patch, 1 if with a shock patch.
%
% sEdgeInd: logical vector, same length as the first dimension of edges. 
% 0 if the edge is of a normal patch, 1 if a shock patch.

function [edges,macroNodes,sNodeInd,sEdgeInd] = genPatchLocs(xBnds, sLocs, pWid, spWid, spMacroDist, nPatch)

%%%% Uncomment the following lines to run as a script
% xBnds = [-pi pi]; %left and right microscale boundaries
% sLocs = [];%[0 :0.5:2.5];%linspace(-pi,pi,11); sLocs = sLocs(2:end-1);%[0]; %(initial) locations of shock(s)
% pWid=0.5; %width of an ordinary patch
% spWid = .5; %width of a shock patch
% spMacroInd = 7; %number of microscale nodes from either edge of a shock patch to 
% %the macroscale node inside it.
% nPatch = 11; %minimum total number of patches (of both types) across domain
% nump = 15; %number of nodes in a normal patch
% nump = 2*floor(nump/2)+1; %current code assumes the core is a single node: make sure nump is odd
% numSp = 15; %number of nodes in a `shock' patch
% numSp = 2*floor(numSp)+1; 
% spH = spWid/(numSp+1); %microscale node distance in shock patch
% pH = pWid/(nump+1); %microscale node distance in normal patch
% spMacroDist = spMacroInd*spH; %distance from either edge of a shock patch to 
% %the macroscale node inside it.
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%

 

sLocs = reshape(sLocs,[],1); %ensure sLocs is a column vector...
sLocs = sort(sLocs); %...in ascending order.

xL = xBnds(1); xR = xBnds(2);

nSP = length(sLocs); %number of shock patches

nP = nPatch-nSP; %number of normal patches

if nP*pWid + nSP*spWid> diff(xBnds)
    error('Requested patches cannot fit in the specified boundaries.')
end
if spWid/2 < spMacroDist
    error('Shock patches are too narrow to accomodate the input spMacroDist, the distance from patch edge to macro-node.')
end

pSpacing = (xR-xL)/(nPatch); %mean spacing between patch edges


%first place `shock patches'
sEdges = [sLocs-spWid/2 sLocs+spWid/2];

edges = [xL xL+pWid; ...
    sEdges; ...
    xR-pWid xR];

if sum( edges(2:end,1)-edges(1:end-1,2)<0 ) > 0
    error('One or more shocks are too close to each other, or to the boundaries; or a patch is located outside the boundaries.')
end


sMacroNodes =  [sEdges(:,1)+spMacroDist sEdges(:,2)-spMacroDist];

macroNodes = sort([xL+pWid/2; sMacroNodes(:); xR-pWid/2]);



remainders = edges(2:end,1)-edges(1:end-1,2); %remaining space(s) to be filled with patches

nRemain = floor(remainders/pSpacing); %number of remaining patches to place... rough estimate.

newNodes = [];

for n=1:length(nRemain)
    nP = nRemain(n);
    lEdge = edges(n,2); %use edges because, unlike nodes, the number of edges is the same in normal&double patches
    rEdge = edges(n+1,1);
    lInd = find(macroNodes - lEdge <0, 1, 'last');
    rInd = find(macroNodes - rEdge >0, 1, 'first');
    lNode = macroNodes(lInd);
    rNode = macroNodes(rInd);
    tmpNodes = linspace(lNode,rNode,nP+2)';
    newNodes = [newNodes; tmpNodes(2:end-1)];
end

macroNodes = sort([macroNodes; newNodes]);
edges = sort([edges; newNodes-pWid/2 newNodes+pWid/2]); %sort goes by column only

sNodeInd=zeros(length(macroNodes),1);
sEdgeInd=zeros(length(edges),1);
for m=1:2*nSP
    sNodeInd(macroNodes == sMacroNodes(m))=1;
end
for m=1:nSP
    sEdgeInd(edges(:,1) == sEdges(m,1))=1;
end





