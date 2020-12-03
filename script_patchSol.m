%Compute patch scheme with one or more double patches.
%Microscale computations done with a holistic (AJR, 2007) discretization
%of Burgers' pde. 
%Trusted solution computed with a numerical quadrature approximation of the
%exact solution to Burgers' equation.
clear; close all

%key parameters
epsilon=0.001; %diffusion coeff in Burgers

%% initial conditions 
iniStr = 'sinusoidal'; %initial condition -sin: model problem M1
% iniStr = 'tanh'; %IC containining shock: model problem M1

%%% sinusoidal IC
[u0,u0int,locPeak,nPatch,sProp] = getIC(iniStr);




xBnds = [-pi pi]; %left and right microscale boundaries
sLocs = [0]; %(initial) locations of shock(s)
pWid=0.01; %width of an ordinary patch
spWid = sProp*pWid; %width of a shock patch
spMacroInd = 2; %number of microscale nodes from either edge of a shock patch to 
%the macroscale node inside it.
lBnd = @(t)[0;0]; %left boundary condition and its derivative
rBnd = @(t)[0;0]; %rbc and its derivative

tSpan = linspace(0, 3, 10);
Gamma = 3; %maximum order of patch coupling - truncated near shocks and boundaries.



nump = 5; %number of nodes in a normal patch
nump = 2*floor(nump/2)+1; %current code assumes the core is a single node: make sure nump is odd
numSp = sProp*(nump+1)-1; %number of nodes in a `shock' patch
% numSp = 2*floor(numSp/2)+1; 
spH = spWid/(numSp+1); %microscale node distance in shock patch
pH = pWid/(nump+1); %microscale node distance in normal patch

spMacroDist = spMacroInd*spH; %distance from either edge of a shock patch to 
%the macroscale node inside it.

%get locations for all patch edges, and the locations of the nodes in each
%patch
[edges,macroNodes,sNodeInd,sEdgeInd] = genPatchLocs(xBnds, sLocs, pWid,spWid, spMacroDist, nPatch);

nPatch = length(edges); %write over (possibly changed) nPatch
nNodes = nPatch+length(sLocs); %two macroscale nodes in each shock patch
%sanity test: nNodes should equal length(macroNodes)


%next learn where (and how many) the independent regions of macroscale nodes are
%The construct intRegions separates patch interpolation across shocks. For example if
%there is one shock, intRegions{1} will hold all nodes to the left of the
%shock, and intRegions{2} will hold all nodes to the right. 
if isempty(sLocs) %no shocks
    numIntRegions = 1;
    intRegions = {(1:nNodes)'};
else %at least one shock
    %shocks separate regions of interpolation
    numIntRegions = 1 + length(sLocs); %how many regions are there?
    intRegions = cell(numIntRegions,1);
    indLocs = find(sNodeInd); %where are the boundary nodes adjacent to a shock?
    intRegions{1} = (1:indLocs(1))';
    intRegions{numIntRegions} = (indLocs(end):nNodes)';
    for m=2:numIntRegions-1
        intRegions{m} = (indLocs(2*m-2):indLocs(2*m-1))';
    end
end

%% get internal pde discretizations for both types of patch
%we are using a holistic discretization - see reference in burgDisc()
gamma=1; %coupling parameter from AJR - always use 1 for simulation
%shock patch discretizations
usp = burgDisc(numSp,spH,gamma,epsilon);
%normal patch discretizations
up = burgDisc(nump,pH,gamma,epsilon);


%(Decide how to) interpolate macroscale values to set patch boundaries
inter = @(xLocs,evalPt) lagInterp(xLocs,evalPt); %lagrangian interpolation


%% Construct the fine-scale grid, and an index of macroscale nodes
xPatch = [];
x2Xind = 0; %initialise at 0, remove after
nxtDist = 0; %distance from current macroscale node to right-hand patch edge
charInd = []; %will be 0s for normal patches, 1s for shock patches
patchInd = cell(nPatch,1); %cell array where each entry is the indices for a particular patch
lastInd = 0;
for n=1:nPatch
    charFun = sEdgeInd(n); %1 if shock patch, 0 otherwise
    tmpN = charFun*numSp + (1-charFun)*nump;
    tmp = linspace(edges(n,1),edges(n,2),tmpN+2)';
    tmp = tmp(2:end-1);
    xPatch = [xPatch; tmp];
    patchInd{n} = lastInd+1:lastInd+length(tmp); lastInd = lastInd+length(tmp);
    if charFun
        x2Xind = [x2Xind; x2Xind(end)+nxtDist+spMacroInd; x2Xind(end)+nxtDist+numSp+1-spMacroInd];
        charInd = [charInd; ones(tmpN,1)];
    else
        x2Xind = [x2Xind; x2Xind(end)+nxtDist+(nump+1)/2];
        charInd = [charInd; zeros(tmpN,1)];
    end
    nxtDist = charFun*(spMacroInd-1) + (1-charFun)*(nump-1)/2;
end
x2Xind = x2Xind(2:end); %remove 0 value

%xPatch is the complete array of patches
%edges holds the locations at which boundary values must be specified
%x2Xind is the same length as nNodes and holds indices of macroscale nodes in xPatch
%sanity test:
% testPatchNodeConsistency = norm(xPatch(x2Xind)-macroNodes) %should give 0

    
%% Recursively construct an ode to simulate in every patch
%Future implementations will be more efficient; but for this version,
%construct each patch ode and their coupling conditions in 
%sequence, then concatenate them. 

%Begin with the left boundary patch
k=1; %current macroscale region
rBndInds = intRegions{k}(abs(intRegions{k}-1) < Gamma+1); %available macroscale nodes
rBndLoc = edges(1,2); %read in boundary locations of the patch
rInter = inter(macroNodes(rBndInds),rBndLoc); %generate a function to interpolate to patch boundary
rMacros = @(u) u(x2Xind(rBndInds)); %get values for macro variables in patch
tmp = @(t,u) up(t,u(patchInd{1}),lBnd(t),[rInter(rMacros(u)); 0]);

%now iterate, adding one patch each time
j=1; %number of nodes passed
for m=2:nPatch-1
    charFun = sEdgeInd(m); %1 if shock patch, 0 otherwise
    
    if charFun
        j=j+1;
        lBndInds = intRegions{k}(abs(intRegions{k}-j) < Gamma+1); %available macroscale nodes
         
        k=k+1; %crossing the shock
         j=j+1;
        rBndInds = intRegions{k}(abs(intRegions{k}-j) < Gamma+1); %available macroscale nodes
        
        lInter = inter(macroNodes(lBndInds),edges(m,1)); %generate a function to interpolate to patch boundary
        rInter = inter(macroNodes(rBndInds),edges(m,2));
        uLMacro = @(u) u(x2Xind(lBndInds));
        uRMacro = @(u) u(x2Xind(rBndInds));
        tmp = @(t,u) [tmp(t,u);
                      usp(t,u(patchInd{m}),[lInter(uLMacro(u));0], [rInter(uRMacro(u));0])];
        
    else
        j=j+1;
        bndInds = intRegions{k}(abs(intRegions{k}-j) < Gamma+1); %available macroscale nodes
        lInter = inter(macroNodes(bndInds),edges(m,1)); %generate a function to interpolate to patch boundary
        rInter = inter(macroNodes(bndInds),edges(m,2));
        uMacro = @(u) u(x2Xind(bndInds));
        tmp = @(t,u) [tmp(t,u);
                      up(t,u(patchInd{m}),[lInter(uMacro(u));0], [rInter(uMacro(u));0])];
    end
end

%now add the right boundary
lBndInds = intRegions{k}(abs(intRegions{k}-nNodes) < Gamma+1); %available macroscale nodes
lBndLoc = edges(end,1); %read in boundary locations of the patch
lInter = inter(macroNodes(lBndInds),lBndLoc); %generate a function to interpolate to patch boundary
lMacros = @(u) u(x2Xind(lBndInds)); %get values for macro variables in patch
%rename the assembled ode to patchOde
patchOde = @(t,u) [tmp(t,u);
              up(t,u(patchInd{nPatch}),[lInter(lMacros(u)); 0],rBnd(t))];


          





%% get full solution
uIC = u0(xPatch);
[t,uSol] = ode15s(patchOde,tSpan,uIC);

[tt,xx] = meshgrid(t,xPatch);

figure(1); hold on
surf(tt,xx,uSol')
xlabel('Time $t$', 'Interpreter', 'Latex')
ylabel('Position $x$', 'Interpreter', 'Latex')
zlabel('Solution $u$', 'Interpreter', 'Latex')
title('Patch Solution')
view(110,35)


%% get trusted solution
tol=5; %sufficient for epsilon<= 0.01
exSol  = trustSol(t,xPatch,epsilon,u0,u0int,locPeak,tol);

figure(2); hold on
surf(tt,xx,exSol')
xlabel('Time $t$', 'Interpreter', 'Latex')
ylabel('Position $x$', 'Interpreter', 'Latex')
zlabel('Solution $u$', 'Interpreter', 'Latex')
title('Trusted Solution')
view(110,35)



pErr = abs(uSol-exSol);
maxErr = max(max(abs(uSol-exSol)))
maxMacroErr = max(max(abs(uSol(:,x2Xind)-exSol(:,x2Xind))))

%% plot errors at specific times
figure(3); 
subplot(2,2,1)
plot(xPatch,uSol(end,:),'o')
hold on
plot(xPatch,exSol(end,:),'o')
xlabel('Position $x$', 'Interpreter', 'Latex')
ylabel('Simulated solution $u$', 'Interpreter', 'Latex')
legend('patch','trusted')
title(['Final time t=' num2str(t(end))])
[iT] = find(max(pErr,[],2) == max(max(pErr)));
subplot(2,2,2)
plot(xPatch,uSol(iT,:),'o')
hold on
plot(xPatch,exSol(iT,:),'o')
legend('patch','trusted')
xlabel('Position $x$', 'Interpreter', 'Latex')
ylabel('Simulated solution $u$', 'Interpreter', 'Latex')
title(['Largest errors at t=' num2str(tSpan(iT)) '.'])
subplot(2,2,3)
semilogy(t, max(pErr,[],2))
xlabel('Time $t$', 'Interpreter', 'Latex')
ylabel('Error (microscale)')
subplot(2,2,4); hold on
semilogy(t, max(pErr(:,x2Xind),[],2))
%semilogy(t, max(pErr(:,x2Xind)1:nPatch,[],2))))
xlabel('Time $t$', 'Interpreter', 'Latex')
ylabel('Error (macroscale)')






    

