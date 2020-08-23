function [FinalX_central_allGps,FinalY_central_allGps,...
    FinalX_UB95_allGps,FinalY_UB95_allGps,FinalX_LB95_allGps,FinalY_LB95_allGps,medium95percError] = ...
    MTcomputation(image2D,NumOfShifts,GridSize,PSFsigma,pixelSize,beta,cRatio,l_KR,gridPtUnit)
%% Feature enhancement
if size(image2D,3) == 3
    image2D = rgb2gray(image2D);
end
if ~isa(image2D,'uint16')
    image2D = im2uint16(image2D./max(image2D(:)));
end
EnhancedProj=imadjust(image2D);
FeatureEnhancedProj=CurLinFeatureEnhancement(double(EnhancedProj),PSFsigma,pixelSize,beta,cRatio);
%% Shifted divisions & Pre-processing B
IMDim=size(FeatureEnhancedProj,1);
Indents=0:round(GridSize/NumOfShifts):(GridSize-1);
IMs_cell=repmat({FeatureEnhancedProj},size(Indents));
GridSize_cell=repmat({GridSize},size(Indents));
Indents_cell=num2cell(Indents);
GaussfiltSigma_cell=repmat({PSFsigma/pixelSize*2},size(Indents));
[AllPt1x,AllPt2x,AllPt1y,AllPt2y,~]=cellfun(@Division2LineSegs,...
    IMs_cell,GridSize_cell,Indents_cell,GaussfiltSigma_cell,'UniformOutput',false);
AllShiftsPt1x=cell2mat(AllPt1x');
AllShiftsPt1y=cell2mat(AllPt1y');
AllShiftsPt2x=cell2mat(AllPt2x');
AllShiftsPt2y=cell2mat(AllPt2y');
%% Grouping
% Remove Edge Noises
dEdge=round(GridSize/NumOfShifts);
xpar=find(AllShiftsPt1x==AllShiftsPt2x & (AllShiftsPt1x<dEdge | AllShiftsPt1x > (IMDim-dEdge)));
ypar=find(AllShiftsPt1y==AllShiftsPt2y & (AllShiftsPt1y<dEdge | AllShiftsPt1y > (IMDim-dEdge)));
EdgeNoise=setdiff(1:length(AllShiftsPt1x),union(xpar,ypar));
AllEndPts=[AllShiftsPt1x AllShiftsPt1y AllShiftsPt2x AllShiftsPt2y];
AllEndPts=AllEndPts(EdgeNoise,:);

% Remove duplicates in original data
[~,ia,~]=unique(AllEndPts,'stable','rows');
Pt1x=AllEndPts(ia,1);
Pt1y=AllEndPts(ia,2);
Pt2x=AllEndPts(ia,3);
Pt2y=AllEndPts(ia,4);

% Calculate Correlation Matrix (1st round)
TotalNum=length(Pt1x);C = combnk(1:TotalNum,2);
Seg1_1x=Pt1x(C(:,1));Seg2_1x=Pt1x(C(:,2));
Seg1_1y=Pt1y(C(:,1));Seg2_1y=Pt1y(C(:,2));
Seg1_2x=Pt2x(C(:,1));Seg2_2x=Pt2x(C(:,2));
Seg1_2y=Pt2y(C(:,1));Seg2_2y=Pt2y(C(:,2));

Angles=atan2(Pt2y-Pt1y,Pt2x-Pt1x);
Angles1=Angles(C(:,1));
Angles2=Angles(C(:,2));

dAngles=wrapTo0halfPi(bsxfun(@minus,Angles2,Angles1));
PolygonArea=arrayfun(@calculatePolygonArea,Seg1_1x,Seg1_1y,Seg1_2x,Seg1_2y,Seg2_1x,Seg2_1y,Seg2_2x,Seg2_2y);
RecArea=sqrt(((Seg1_2x-Seg1_1x).^2+(Seg1_2y-Seg1_1y).^2).*((Seg2_2x-Seg2_1x).^2+(Seg2_2y-Seg2_1y).^2));
PolygonAreaRatio=PolygonArea./RecArea;

D11p=sqrt(((Seg2_1x-Seg1_1x).^2+(Seg2_1y-Seg1_1y).^2));
D12p=sqrt(((Seg2_2x-Seg1_1x).^2+(Seg2_2y-Seg1_1y).^2));
D21p=sqrt(((Seg1_2x-Seg2_1x).^2+(Seg1_2y-Seg2_1y).^2));
D22p=sqrt(((Seg1_2x-Seg2_2x).^2+(Seg1_2y-Seg2_2y).^2));
minCrossDist=min([D11p D12p D21p D22p],[],2);
minCrossDist(minCrossDist==0)=eps;% to avoid 0
Len1=sqrt(((Seg1_2x-Seg1_1x).^2+(Seg1_2y-Seg1_1y).^2));
Len2=sqrt(((Seg2_2x-Seg2_1x).^2+(Seg2_2y-Seg2_1y).^2));
avgLen=mean([Len1 Len2],2);
LenRatio=minCrossDist./avgLen;

nei=find(LenRatio<1);
ang=find(dAngles<0.25*pi);
are=find(PolygonAreaRatio<0.5);
sel=intersect(intersect(ang,are),nei);

k=LenRatio.*PolygonAreaRatio.*1./(cos(dAngles).^2);
A=exp(-k);
mu=mean(A(sel));
sigma=std(A(sel),1);
h=-mean(A(sel).*log(A(sel)));
ep=mu-sigma+h;
Asel=find(A>ep);
FirstRound=intersect(sel,Asel);

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
revisedGpNum=(min(uniqueGps):max(uniqueGps))';
groupInd=revisedGpNum(ic);
%% Calculate averaged curves and errors for each segment group
FinalX_central_allGps=cell(max(groupInd),1);
FinalY_central_allGps=cell(max(groupInd),1);
FinalX_UB95_allGps=cell(max(groupInd),1);
FinalY_UB95_allGps=cell(max(groupInd),1);
FinalX_LB95_allGps=cell(max(groupInd),1);
FinalY_LB95_allGps=cell(max(groupInd),1);
medium95percError=zeros(max(groupInd),1);

for i=1:max(groupInd)
    % 0. extract a group
    CurrGroupInd=find(groupInd==i);
    P1x=Pt1x(CurrGroupInd);P2x=Pt2x(CurrGroupInd);
    P1y=Pt1y(CurrGroupInd);P2y=Pt2y(CurrGroupInd);
    
    if length(P1x) < 3
        continue;
    end
    
    % 1. get the principal direction
    EndPts=[P1x P1y;P2x P2y];
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
    X=[ptsX(:) ptsY(:)];
    X_cell=mat2cell(X,repmat(length(P1x),1,n),2);
    
    principalCurve_cell=cellfun(@(x) computePC(x, l_KR, gridPtUnit),X_cell,'UniformOutput',false);
    principalCurve_cell_T=cellfun(@(x) coeffEndPts\(x-mean(EndPts))', principalCurve_cell, 'UniformOutput',false);% transform coordinates
    
    % 3. compute Curve Average with extrapolation
    [FinalX_T, FinalY_central_T, FinalY_95up_T, FinalY_95low_T]=computeCurveAverage(principalCurve_cell_T,xmin_T,xmax_T,gridPtUnit);
    
    error = median((FinalY_95up_T-FinalY_95low_T)./2);
    
    if error > 1
        continue; % remove curves with very large error
    end
    
    % 4. transform back to standard coordinates
    FinalX_centralandFinalY_central=coeffEndPts*[FinalX_T;FinalY_central_T]+mean(EndPts)';
    FinalX_UBandFinalY_UB95=coeffEndPts*[FinalX_T;FinalY_95up_T]+mean(EndPts)';
    FinalX_LBandFinalY_LB95=coeffEndPts*[FinalX_T;FinalY_95low_T]+mean(EndPts)';
    
    FinalX_central=FinalX_centralandFinalY_central(1,:);
    FinalY_central=FinalX_centralandFinalY_central(2,:);
    FinalX_UB95=FinalX_UBandFinalY_UB95(1,:);
    FinalY_UB95=FinalX_UBandFinalY_UB95(2,:);
    FinalX_LB95=FinalX_LBandFinalY_LB95(1,:);
    FinalY_LB95=FinalX_LBandFinalY_LB95(2,:);
    
    FinalX_central_allGps(i)={FinalX_central};
    FinalY_central_allGps(i)={FinalY_central};
    FinalX_UB95_allGps(i)={FinalX_UB95};
    FinalY_UB95_allGps(i)={FinalY_UB95};
    FinalX_LB95_allGps(i)={FinalX_LB95};
    FinalY_LB95_allGps(i)={FinalY_LB95};
    medium95percError(i)=error;
end

FinalX_central_allGps=FinalX_central_allGps(~cellfun('isempty',FinalX_central_allGps));
FinalY_central_allGps=FinalY_central_allGps(~cellfun('isempty',FinalY_central_allGps));
FinalX_UB95_allGps=FinalX_UB95_allGps(~cellfun('isempty',FinalX_UB95_allGps));
FinalY_UB95_allGps=FinalY_UB95_allGps(~cellfun('isempty',FinalY_UB95_allGps));
FinalX_LB95_allGps=FinalX_LB95_allGps(~cellfun('isempty',FinalX_LB95_allGps));
FinalY_LB95_allGps=FinalY_LB95_allGps(~cellfun('isempty',FinalY_LB95_allGps));
medium95percError=medium95percError(medium95percError~=0);

