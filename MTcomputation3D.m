function [FinalX_central_allGps,FinalY_central_allGps,FinalZ_central_allGps,...
    FinalX_B95_allGps, FinalY_B95_allGps, FinalZ_B95_allGps,medium95percError] = ...
    MTcomputation3D(image3D,NumOfShifts,GridSize,PSFsigma,pixelSize,alpha,beta,cRatio,l_KR,gridPtUnit,num_seg_thresh)
%% Feature enhancement requires a grayscale input

if ~isa(image3D,'uint16')
    image3D = im2uint16(image3D./max(image3D(:)));
end
EnhancedIM=imadjustn(image3D);
FeatureEnhancedIM=CurLinFeatureEnhancement3D_ellip(double(EnhancedIM),PSFsigma,pixelSize,alpha,beta,cRatio);

disp('Image enhanced');
%% Shifted divisions & Pre-processing B (rectangular)
if isscalar(PSFsigma)
    PSFsigma = [PSFsigma PSFsigma PSFsigma];
elseif length(PSFsigma) == 2
    PSFsigma = [PSFsigma(1) PSFsigma(1) PSFsigma(2)];
end
if isscalar(pixelSize)
    pixelSize = [pixelSize pixelSize pixelSize];
elseif length(pixelSize)==2
    pixelSize = [pixelSize(1) pixelSize(1) pixelSize(2)];
end
if isscalar(GridSize)
    GridSize = [GridSize GridSize size(image3D,3)];
elseif length(GridSize) == 2
    GridSize = [GridSize(1) GridSize(1) GridSize(2)]; 
end
if isscalar(NumOfShifts)
    NumOfShifts = [NumOfShifts NumOfShifts 1];
elseif length(NumOfShifts) == 2
    NumOfShifts = [NumOfShifts(1) NumOfShifts(1) NumOfShifts(2)]; 
end

IMDim=size(image3D);

if (GridSize(3) == IMDim(3)) || (NumOfShifts(3) == 1)
    [ia,ib] = meshgrid(0:round(GridSize(1)/NumOfShifts(1)):(GridSize(1)-1),...
    0:round(GridSize(2)/NumOfShifts(2)):(GridSize(2)-1));
    ic = zeros(size(ia));
else
[ia,ib,ic] = meshgrid(0:round(GridSize(1)/NumOfShifts(1)):(GridSize(1)-1),...
    0:round(GridSize(2)/NumOfShifts(2)):(GridSize(2)-1),...
    0:round(GridSize(3)/NumOfShifts(3)):(GridSize(3)-1));
end
Indents = [ia(:) ib(:) ic(:)]; 

Indents = Indents(Indents(:,1)==Indents(:,2),:);
scale2PSFratio = 1/sqrt(2);

%%%%%%%%%%%%% use cellfun %%%%%%%%%%%%%%%%%%%%%%%
IMs_cell=repmat({FeatureEnhancedIM},size(Indents,1),1);
GridSize_cell=repmat({GridSize},size(Indents,1),1);
Indents_cell=mat2cell(Indents,ones(size(Indents,1),1),3); 
GaussfiltSigma_cell=repmat({PSFsigma./pixelSize*scale2PSFratio},size(Indents,1),1);
[AllPt1x,AllPt2x,AllPt1y,AllPt2y,AllPt1z,AllPt2z]=cellfun(@Division2LineSegs3D,...
    IMs_cell,GridSize_cell,Indents_cell,GaussfiltSigma_cell,'UniformOutput',false);
%%%%%%%%%%%%% use cellfun %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% rewrite with a for loop %%%%%%%%%%%%%%%%%%
% AllPt1x = cell(size(Indents,1),1);
% AllPt2x = cell(size(Indents,1),1);
% AllPt1y = cell(size(Indents,1),1);
% AllPt2y = cell(size(Indents,1),1);
% AllPt1z = cell(size(Indents,1),1);
% AllPt2z = cell(size(Indents,1),1);
% 
% for dl = 1:size(Indents,1)
%     [tAllPt1x,tAllPt2x,tAllPt1y,tAllPt2y,tAllPt1z,tAllPt2z]=Division2LineSegs3D(...
%     FeatureEnhancedIM,GridSize,Indents(dl,:),PSFsigma./pixelSize*scale2PSFratio);
%     AllPt1x(dl) = {tAllPt1x};
%     AllPt2x(dl) = {tAllPt2x};
%     AllPt1y(dl) = {tAllPt1y};
%     AllPt2y(dl) = {tAllPt2y};
%     AllPt1z(dl) = {tAllPt1z};
%     AllPt2z(dl) = {tAllPt2z};
%     dl
% end
%%%%%%%%%%%% rewrite with a for loop %%%%%%%%%%%%%%%%%%


AllShiftsPt1x=cell2mat(AllPt1x(:));
AllShiftsPt1y=cell2mat(AllPt1y(:));
AllShiftsPt1z=cell2mat(AllPt1z(:));
AllShiftsPt2x=cell2mat(AllPt2x(:));
AllShiftsPt2y=cell2mat(AllPt2y(:));
AllShiftsPt2z=cell2mat(AllPt2z(:));

if isempty(AllShiftsPt1x)
    FinalX_central_allGps=[];
    FinalY_central_allGps=[];
    FinalZ_central_allGps=[];
    FinalX_B95_allGps=[];
    FinalY_B95_allGps=[];
    FinalZ_B95_allGps=[];
    medium95percError=[];
    disp('No feature detected');
    return;
end

disp('Line segments detected');
%% Grouping
% Remove Edge Noise
dEdge=GridSize./NumOfShifts;
xpar=find(AllShiftsPt1y==AllShiftsPt2y & AllShiftsPt1z==AllShiftsPt2z & ...
    (AllShiftsPt1y<dEdge(2) | AllShiftsPt1y > (IMDim(2)-dEdge(2))) & ...
    (AllShiftsPt1z<dEdge(3) | AllShiftsPt1z > (IMDim(3)-dEdge(3))));
ypar=find(AllShiftsPt1x==AllShiftsPt2x & AllShiftsPt1z==AllShiftsPt2z & ...
    (AllShiftsPt1x<dEdge(1) | AllShiftsPt1x > (IMDim(1)-dEdge(1))) & ...
    (AllShiftsPt1z<dEdge(3) | AllShiftsPt1z > (IMDim(3)-dEdge(3))));
zpar=find(AllShiftsPt1x==AllShiftsPt2x & AllShiftsPt1y==AllShiftsPt2y & ...
    (AllShiftsPt1x<dEdge(1) | AllShiftsPt1x > (IMDim(1)-dEdge(1))) & ...
    (AllShiftsPt1y<dEdge(2) | AllShiftsPt1y > (IMDim(2)-dEdge(2))));

EdgeNoise=setdiff(1:length(AllShiftsPt1x),union(union(xpar,ypar),zpar));
AllEndPts=[AllShiftsPt1x AllShiftsPt1y AllShiftsPt1z AllShiftsPt2x AllShiftsPt2y AllShiftsPt2z];
AllEndPts=AllEndPts(EdgeNoise,:);

% Remove duplicates in original data
[~,ia,~]=unique(AllEndPts,'stable','rows');
if size(ia,1) < 3
    FinalX_central_allGps=[];
    FinalY_central_allGps=[];
    FinalZ_central_allGps=[];
    FinalX_B95_allGps=[];
    FinalY_B95_allGps=[];
    FinalZ_B95_allGps=[];
    medium95percError=[];
    disp('No feature detected');
    return;
end
Pt1 = AllEndPts(ia,1:3);
Pt2 = AllEndPts(ia,4:6);

% Calculate Correlation Matrix (1st round)
TotalNum=size(Pt1,1);C = combnk(1:TotalNum,2);
Seg1_P1 = Pt1(C(:,1),:);Seg2_P1 = Pt1(C(:,2),:);
Seg1_P2 = Pt2(C(:,1),:);Seg2_P2 = Pt2(C(:,2),:);
Seg1_len = vecnorm((Seg1_P1-Seg1_P2)')';
Seg2_len = vecnorm((Seg2_P1-Seg2_P2)')';
avg_len = (Seg1_len+Seg2_len)./2;

% initialize final parameters
vertd = zeros(size(Seg1_P1,1),1);
interAngle = zeros(size(Seg1_P1,1),1);
LenRatio=zeros(size(Seg1_P1,1),1);
PolygonAreaRatio = zeros(size(Seg1_P1,1),1);

% First decide if coplanar => the fourth parameter (vertd) should be set to 0
AminusD = Seg1_P1 - Seg2_P2;
BminusD = Seg1_P2 - Seg2_P2;
CminusD = Seg2_P1 - Seg2_P2;
V = abs(dot(AminusD,cross(BminusD,CminusD,2),2))/6;
coplanarInd = (V==0);

% Calculate the moving vector only for V~=0 case
% p1=r1+t1e1; p2=r2+t2e2; n=e1 crossX e2; d = n(r1-r2)/|n|; 
e_Seg1 = (Seg1_P1-Seg1_P2) ./ Seg1_len;
e_Seg2 = (Seg2_P1-Seg2_P2) ./ Seg2_len;
vmove = cross(e_Seg1, e_Seg2,2);
dpz = dot(vmove,(Seg1_P1-Seg2_P1),2)./(vecnorm(vmove')');
dpz(coplanarInd) = 0;
noncopla = ~coplanarInd;
vertd(noncopla) = abs(dpz(noncopla))./(avg_len(noncopla));

% Move Seg2
Seg2_P1_moved = Seg2_P1;
Seg2_P2_moved = Seg2_P2;
Seg2_P1_moved(noncopla,:) = Seg2_P1(noncopla,:)+(vmove(noncopla,:)./...
    (vecnorm(vmove(noncopla,:)')').*dpz(noncopla));
Seg2_P2_moved(noncopla,:) = Seg2_P2(noncopla,:)+(vmove(noncopla,:)./...
    (vecnorm(vmove(noncopla,:)')').*dpz(noncopla));

% Transfer axis to the same plane for the non-coplanar ones
localOrigin = Seg1_P1';
localXvector = Seg1_P1-Seg1_P2;
n1 = cross(Seg1_P1-Seg1_P2, Seg2_P1_moved-Seg2_P2_moved,2); n1_norm = vecnorm(n1')';
n2 = cross(Seg1_P1-Seg1_P2, Seg1_P1-Seg2_P1_moved,2); n2_norm = vecnorm(n2')';
n3 = cross(Seg1_P1-Seg1_P2, Seg1_P1-Seg2_P2_moved,2); n3_norm = vecnorm(n3')';
[maxnorm, n1_max] = max([n1_norm, n2_norm, n3_norm],[],2);
colinearInd = (maxnorm==0); 
localZvector = zeros(size(localXvector));
localZvector(n1_max==1,:) = n1(n1_max==1,:);
localZvector(n1_max==2,:) = n2(n1_max==2,:);
localZvector(n1_max==3,:) = n3(n1_max==3,:);

localYvector = cross(localXvector,localZvector,2);
localAxes_x = localXvector./(vecnorm(localXvector')');
localAxes_y = localYvector./(vecnorm(localYvector')');
localAxes_z = localZvector./(vecnorm(localZvector')');
localAxes = cat(2,reshape(localAxes_x',3,1,size(localAxes_x,1)),...
    reshape(localAxes_y',3,1,size(localAxes_y,1)),...
    reshape(localAxes_z',3,1,size(localAxes_z,1))); 

Seg1_P1_transf= global2localcoord(Seg1_P1(~colinearInd,:)',...
    'rr',localOrigin(:,~colinearInd),localAxes(:,:,~colinearInd))';
Seg1_P2_transf = global2localcoord(Seg1_P2(~colinearInd,:)',...
    'rr',localOrigin(:,~colinearInd),localAxes(:,:,~colinearInd))';
Seg2_P1_transf = global2localcoord(Seg2_P1_moved(~colinearInd,:)',...
    'rr',localOrigin(:,~colinearInd),localAxes(:,:,~colinearInd))';
Seg2_P2_transf = global2localcoord(Seg2_P2_moved(~colinearInd,:)',...
    'rr',localOrigin(:,~colinearInd),localAxes(:,:,~colinearInd))';
% indexing of Seg1_P1_transf is different from Seg1_P1 and LenRatio, etc.

% Calculate the four parameters: angle
Seg1_P1_rSpace = [(Seg1_P1_transf(:,1)-0.5).*pixelSize(1) ...
    (Seg1_P1_transf(:,2)-0.5).*pixelSize(2) (Seg1_P1_transf(:,3)-0.5).*pixelSize(3)];
Seg1_P2_rSpace = [(Seg1_P2_transf(:,1)-0.5).*pixelSize(1) ...
    (Seg1_P2_transf(:,2)-0.5).*pixelSize(2) (Seg1_P2_transf(:,3)-0.5).*pixelSize(3)];
Seg2_P1_rSpace = [(Seg2_P1_transf(:,1)-0.5).*pixelSize(1) ...
    (Seg2_P1_transf(:,2)-0.5).*pixelSize(2) (Seg2_P1_transf(:,3)-0.5).*pixelSize(3)];
Seg2_P2_rSpace = [(Seg2_P2_transf(:,1)-0.5).*pixelSize(1) ...
    (Seg2_P2_transf(:,2)-0.5).*pixelSize(2) (Seg2_P2_transf(:,3)-0.5).*pixelSize(3)];

cosAngle = dot((Seg1_P1_rSpace-Seg1_P2_rSpace),(Seg2_P1_rSpace-Seg2_P2_rSpace),...
    2)./(vecnorm((Seg1_P1_rSpace-Seg1_P2_rSpace)')')./(vecnorm((Seg2_P1_rSpace-Seg2_P2_rSpace)')');
cosAngle(cosAngle<-1) = -1; % adjust value to make sure no complex results
cosAngle(cosAngle>1) = 1;
dAngles=wrapTo0halfPi(acos(cosAngle));
interAngle(~colinearInd) = dAngles;

% Calculate the four parameters: distance
D11p=vecnorm((Seg1_P1_transf-Seg2_P1_transf)')';
D12p=vecnorm((Seg1_P1_transf-Seg2_P2_transf)')';
D21p=vecnorm((Seg1_P2_transf-Seg2_P1_transf)')';
D22p=vecnorm((Seg1_P2_transf-Seg2_P2_transf)')';
minCrossDist=min([D11p D12p D21p D22p],[],2);
minCrossDist(minCrossDist==0)=eps;% to avoid 0
Len1=vecnorm((Seg1_P2_transf-Seg1_P1_transf)')';
Len2=vecnorm((Seg2_P2_transf-Seg2_P1_transf)')';
avgLen=mean([Len1 Len2],2);
LenRatio(~colinearInd)=minCrossDist./avgLen;

% added 20200610 need to calculate lenRatio for colinear pairs
D11p_cl=vecnorm((Seg1_P1(colinearInd,:)-Seg2_P1(colinearInd,:))')';
D12p_cl=vecnorm((Seg1_P1(colinearInd,:)-Seg2_P2(colinearInd,:))')';
D21p_cl=vecnorm((Seg1_P2(colinearInd,:)-Seg2_P1(colinearInd,:))')';
D22p_cl=vecnorm((Seg1_P2(colinearInd,:)-Seg2_P2(colinearInd,:))')';
minCrossDist_cl=min([D11p_cl D12p_cl D21p_cl D22p_cl],[],2);
minCrossDist_cl(minCrossDist_cl==0)=eps;% to avoid 0
Len1_cl=vecnorm((Seg1_P2(colinearInd,:)-Seg1_P1(colinearInd,:))')';
Len2_cl=vecnorm((Seg2_P2(colinearInd,:)-Seg2_P1(colinearInd,:))')';
avgLen_cl=mean([Len1_cl Len2_cl],2);
LenRatio(colinearInd) = minCrossDist_cl./avgLen_cl;

% Calculate the four parameters: area
PolygonArea=arrayfun(@calculatePolygonArea,Seg1_P1_transf(:,1),Seg1_P1_transf(:,2),...
    Seg1_P2_transf(:,1),Seg1_P2_transf(:,2),Seg2_P1_transf(:,1),Seg2_P1_transf(:,2),...
    Seg2_P2_transf(:,1),Seg2_P2_transf(:,2));
RecArea=vecnorm((Seg1_P2_transf-Seg1_P1_transf)').*vecnorm((Seg2_P2_transf-Seg2_P1_transf)');
PolygonAreaRatio(~colinearInd)=PolygonArea./RecArea(:);

% Put everything together (back with the original number of pairs)
nei=find(LenRatio<1);
ang=find(interAngle<pi/4);
are=find(PolygonAreaRatio<0.5);
ver=find(vertd<0.5);
sel=intersect(intersect(intersect(ang,are),nei),ver);

k=LenRatio.*PolygonAreaRatio.*1./(cos(interAngle).^2).*(vertd+1);
A=exp(-k);
mu=mean(A(sel));
sigma=std(A(sel),1);
h=-mean(A(sel).*log(A(sel)));
ep=mu-sigma+h;
Asel=find(A>=ep);
FirstRound=intersect(sel,Asel);

if isempty(FirstRound)
    FinalX_central_allGps=[];
    FinalY_central_allGps=[];
    FinalZ_central_allGps=[];
    FinalX_B95_allGps=[];
    FinalY_B95_allGps=[];
    FinalZ_B95_allGps=[];
    medium95percError=[];
    disp('No feature detected');
    return;
end  

% Group line segment pairs sharing a common segment
groupInd=zeros(TotalNum,1);j=1;
for i=1:length(FirstRound)
    if groupInd(C(FirstRound(i),1)) == 0 && groupInd(C(FirstRound(i),2)) == 0
        groupInd(C(FirstRound(i),1))=j; % both belong to the j-th group
        groupInd(C(FirstRound(i),2))=j;j=j+1;
    elseif groupInd(C(FirstRound(i),1)) == 0 && groupInd(C(FirstRound(i),2)) ~= 0
        groupInd(C(FirstRound(i),1))=groupInd(C(FirstRound(i),2));
    elseif groupInd(C(FirstRound(i),1)) ~= 0 && groupInd(C(FirstRound(i),2)) == 0
        groupInd(C(FirstRound(i),2))=groupInd(C(FirstRound(i),1));
    elseif groupInd(C(FirstRound(i),1))~=groupInd(C(FirstRound(i),2))
        ind=groupInd==groupInd(C(FirstRound(i),2));
        groupInd(ind)=groupInd(C(FirstRound(i),1));
    end
end
[uniqueGps,~,ic]=unique(groupInd);
if ismember(0,uniqueGps)
    revisedGpNum=(uniqueGps(1):uniqueGps(end))';
else
    revisedGpNum=(uniqueGps(1):uniqueGps(end))'-uniqueGps(1)+1;
end
groupInd=revisedGpNum(ic);

disp('Line segments grouped');

%% Calculate averaged curves and errors for each segment group
FinalX_central_allGps=cell(max(groupInd),1);
FinalY_central_allGps=cell(max(groupInd),1);
FinalZ_central_allGps=cell(max(groupInd),1);
FinalX_B95_allGps=cell(max(groupInd),1);
FinalY_B95_allGps=cell(max(groupInd),1);
FinalZ_B95_allGps=cell(max(groupInd),1);
medium95percError=zeros(max(groupInd),1);

for i=1:max(groupInd)
    % 0. extract a group
    CurrGroupInd=find(groupInd==i);
    P1x=Pt1(CurrGroupInd,1);P2x=Pt2(CurrGroupInd,1);
    P1y=Pt1(CurrGroupInd,2);P2y=Pt2(CurrGroupInd,2);
    P1z=Pt1(CurrGroupInd,3);P2z=Pt2(CurrGroupInd,3);
       
    if length(P1x) < num_seg_thresh
        continue; % short feature filtering;
    end
    
    % 1. get the principal line direction using svd
    EndPts=[P1x P1y P1z;P2x P2y P2z];
    EndPts_centered = EndPts-mean(EndPts);
    coeffEndPts = pca(EndPts_centered);
    EndPts_T=coeffEndPts\(EndPts_centered'); 
    xmin_T=min(EndPts_T(1,:));
    xmax_T=max(EndPts_T(1,:));
    
    % 2. use cellfun to get many of these curves
    n=200;
    lambda=rand(length(P1x)*n,1);
    ptsX=lambda.*repmat(P1x,n,1)+(1-lambda).*repmat(P2x,n,1);
    ptsY=lambda.*repmat(P1y,n,1)+(1-lambda).*repmat(P2y,n,1);
    ptsZ=lambda.*repmat(P1z,n,1)+(1-lambda).*repmat(P2z,n,1);
    X=[ptsX(:) ptsY(:) ptsZ(:)];
    X_cell=mat2cell(X,repmat(length(P1x),1,n),3);
    
    principalCurve_cell=cellfun(@(x) computePC_3D(x, l_KR, gridPtUnit),X_cell,'UniformOutput',false);
    principalCurve_cell_T=cellfun(@(x) coeffEndPts\(x-mean(EndPts))', principalCurve_cell, 'UniformOutput',false);% transform coordinates
    
    % 3. compute Curve Average with extrapolation
    [FinalX_T, FinalY_central_T, FinalY_B95_T, FinalZ_central_T, ...
        FinalZ_B95_T]=computeCurveAverage_3D_roundU(principalCurve_cell_T,...
        xmin_T,xmax_T,gridPtUnit);
         
    % 4. transform back to standard coordinates
    FinalX_centralandFinalYZ_central=coeffEndPts*[FinalX_T;FinalY_central_T;FinalZ_central_T]+mean(EndPts)';
    FinalX_central=FinalX_centralandFinalYZ_central(1,:);
    FinalY_central=FinalX_centralandFinalYZ_central(2,:);
    FinalZ_central=FinalX_centralandFinalYZ_central(3,:);
    
    FinalX_T_rep = repmat(FinalX_T',1,size(FinalY_B95_T,2));
    FinalX_andFinalYZ95=coeffEndPts*[FinalX_T_rep(:)';FinalY_B95_T(:)';FinalZ_B95_T(:)']+mean(EndPts)';
    FinalX_B95=reshape(FinalX_andFinalYZ95(1,:),size(FinalY_B95_T,1),size(FinalY_B95_T,2));
    FinalY_B95=reshape(FinalX_andFinalYZ95(2,:),size(FinalY_B95_T,1),size(FinalY_B95_T,2));
    FinalZ_B95=reshape(FinalX_andFinalYZ95(3,:),size(FinalY_B95_T,1),size(FinalY_B95_T,2));
    
    error = median(mean(sqrt((FinalY_B95(:,1:size(FinalY_B95,2)/2)-FinalY_B95(:,size(FinalY_B95,2)/2+1:end)).^2+...
         (FinalZ_B95(:,1:size(FinalZ_B95,2)/2)-FinalZ_B95(:,size(FinalZ_B95,2)/2+1:end)).^2)./2,2));

    if error > 2
        continue; % Large Error filtering 
    end

    FinalX_central_allGps(i)={FinalX_central};
    FinalY_central_allGps(i)={FinalY_central};
    FinalZ_central_allGps(i)={FinalZ_central};
    FinalX_B95_allGps(i)={FinalX_B95};
    FinalY_B95_allGps(i)={FinalY_B95};
    FinalZ_B95_allGps(i)={FinalZ_B95};
    medium95percError(i)=error;
end

FinalX_central_allGps=FinalX_central_allGps(~cellfun('isempty',FinalX_central_allGps));
FinalY_central_allGps=FinalY_central_allGps(~cellfun('isempty',FinalY_central_allGps));
FinalZ_central_allGps=FinalZ_central_allGps(~cellfun('isempty',FinalZ_central_allGps));
FinalX_B95_allGps=FinalX_B95_allGps(~cellfun('isempty',FinalX_B95_allGps));
FinalY_B95_allGps=FinalY_B95_allGps(~cellfun('isempty',FinalY_B95_allGps));
FinalZ_B95_allGps=FinalZ_B95_allGps(~cellfun('isempty',FinalZ_B95_allGps));
medium95percError=medium95percError(medium95percError~=0);