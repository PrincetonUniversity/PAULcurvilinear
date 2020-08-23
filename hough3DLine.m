function [H,Hflat,az_bin,el_bin,xp_bin,yp_bin,Lines] = hough3DLine(BW_or_PtCloud, Threshold_ratio, xp_yp_binsize, minLength, numpeaks, FillGap)
% function similar to the combination of hough, houghpeaks, houghlines
% select the maximum element in HT matrix, no maxPeaks limit
% output structs called lines

if nargin < 2
    Threshold_ratio = 0.5; xp_yp_binsize = 1; minLength = 7; numpeaks = 3; FillGap=4;
elseif nargin < 3
    xp_yp_binsize = 1; minLength = 7; numpeaks = 3; FillGap=4;
elseif nargin < 4
    minLength = 7; numpeaks = 3; FillGap=4;
elseif nargin < 5
    numpeaks = 3; FillGap=4;
elseif nargin < 6
    FillGap=4;
end

if islogical(BW_or_PtCloud)
    [I, J, K] = ind2sub(size(BW_or_PtCloud),find(BW_or_PtCloud));
    PtCloud = [I J K];
elseif isnumeric(BW_or_PtCloud)
    PtCloud = BW_or_PtCloud;
end
if size(PtCloud,1) <= 1
   H=[];Hflat=[];az_bin=[];el_bin=[];xp_bin=[];yp_bin=[];Lines=[];
   return;
end

% get all clusters; if two are separates by less than FillGap, join them
[idx1,idx2,idx3] = ind2sub(size(BW_or_PtCloud),1:length(BW_or_PtCloud(:)));
IDX = [idx1' idx2' idx3'];
CC = bwconncomp(BW_or_PtCloud);
ccpair = combnk(1:length(CC.PixelIdxList),2);
joinPairIdx = zeros(size(ccpair,1),1);
for q = 1:size(ccpair,1)
    D = pdist2(IDX(CC.PixelIdxList{ccpair(q,1)},:),...
        IDX(CC.PixelIdxList{ccpair(q,2)},:));
    if min(D(:)) <= FillGap
        joinPairIdx(q) = 1;
    end
end
posPair = find(joinPairIdx);

if ~isempty(posPair)
    groupInd=zeros(length(CC.PixelIdxList),1);j=1;
    for i=1:length(posPair)
        if groupInd(ccpair(posPair(i),1)) == 0 && groupInd(ccpair(posPair(i),2)) == 0
            groupInd(ccpair(posPair(i),1))=j; % both belong to the j-th group
            groupInd(ccpair(posPair(i),2))=j;j=j+1;
        elseif groupInd(ccpair(posPair(i),1)) == 0 && groupInd(ccpair(posPair(i),2)) ~= 0
            groupInd(ccpair(posPair(i),1))=groupInd(ccpair(posPair(i),2));
        elseif groupInd(ccpair(posPair(i),1)) ~= 0 && groupInd(ccpair(posPair(i),2)) == 0
            groupInd(ccpair(posPair(i),2))=groupInd(ccpair(posPair(i),1));
        elseif groupInd(ccpair(posPair(i),1))~=groupInd(ccpair(posPair(i),2))
            ind=groupInd==groupInd(ccpair(posPair(i),2));
            groupInd(ind)=groupInd(ccpair(posPair(i),1));
        end
    end
    [uniqueGps,~,ic]=unique(groupInd);
    if ismember(0,uniqueGps)
        revisedGpNum=(uniqueGps(1):uniqueGps(end))';
    else
        revisedGpNum=(uniqueGps(1):uniqueGps(end))'-uniqueGps(1)+1;
    end
    groupInd=revisedGpNum(ic);
    PtCloudGroup={};
    singleElementInd = find(groupInd==0);
    PtCloudGroup(1:length(singleElementInd)) = CC.PixelIdxList(singleElementInd);
    for g = 1:max(groupInd)
        CurrGroupInd=groupInd==i;
        currSegGroup = cell2mat(CC.PixelIdxList(CurrGroupInd)');
        PtCloudGroup = [PtCloudGroup {currSegGroup}];
    end
else
    PtCloudGroup = CC.PixelIdxList;
end

% Do HT3D on each cluster of voxels
lines_all = struct('point1',[],'point2',[],'az',[],'el',[],'xp',[],'yp',[],'H',[]);
for v = 1:length(PtCloudGroup)
    
    [currI, currJ, currK] = ind2sub(size(BW_or_PtCloud),PtCloudGroup{v});
    currPtCloud = [currI currJ currK];

% get all point pairs
indcomb = combnk(1:size(currPtCloud,1),2);
Pt1 = currPtCloud(indcomb(:,1),:);
Pt2 = currPtCloud(indcomb(:,2),:);
az = zeros(size(Pt1,1),1);
el = zeros(size(Pt1,1),1);
xp = zeros(size(Pt1,1),1);
yp = zeros(size(Pt1,1),1);
for i=1:size(Pt1,1)
    [az(i), el(i), xp(i), yp(i)] = convertTo4ParmLineRepDeg_cart(Pt1(i,:),Pt2(i,:));
end
Hflat = [az el xp yp];

if size(Hflat,1) == 1
    H=[];az_bin=[];el_bin=[];xp_bin=[];yp_bin=[];Lines=[];
    if DistanceBtPoints(currPtCloud(1,:),currPtCloud(2,:)) >= minLength
        Lines = struct('point1',currPtCloud(1,:),'point2',currPtCloud(2,:),'az',az,...
        'el',el,'xp',xp,'yp',yp,'H',1);
    end
    return;
end

% binning
az_bin = (floor(min(az))-0.5):1:(ceil(max(az))+0.5);if isscalar(az_bin); az_bin = [round(az_bin)-0.5 round(az_bin)+0.5]; end
el_bin = (floor(min(el))-0.5):1:(ceil(max(el))+0.5);if isscalar(el_bin); el_bin = [round(el_bin)-0.5 round(el_bin)+0.5]; end
xp_bin = min(xp):xp_yp_binsize:(max(xp)+xp_yp_binsize);if isscalar(xp_bin); xp_bin = [xp_bin-xp_yp_binsize/2 xp_bin+xp_yp_binsize/2]; end
yp_bin = min(yp):xp_yp_binsize:(max(yp)+xp_yp_binsize);if isscalar(yp_bin); yp_bin = [yp_bin-xp_yp_binsize/2 yp_bin+xp_yp_binsize/2]; end

H = zeros(length(az_bin)-1,length(el_bin)-1,length(xp_bin)-1,length(yp_bin)-1,'uint16');
H_ptPairInd = cell(length(az_bin)-1,length(el_bin)-1,length(xp_bin)-1,length(yp_bin)-1);
for i=1:size(Hflat,1)
    dim1=find(histcounts(Hflat(i,1),az_bin)');
    dim2=find(histcounts(Hflat(i,2),el_bin)');
    dim3=find(histcounts(Hflat(i,3),xp_bin)');
    dim4=find(histcounts(Hflat(i,4),yp_bin)');
    H(dim1,dim2,dim3,dim4)=H(dim1,dim2,dim3,dim4)+1;
    H_ptPairInd(dim1,dim2,dim3,dim4) = {uint16([i H_ptPairInd{dim1,dim2,dim3,dim4}])};
end

if max(H(:)) == 0
    H=[];az_bin=[];el_bin=[];xp_bin=[];yp_bin=[];Lines = [];
    return;
end

% selecting peaks with minimim neighboring distance
Htemp = H;
H_Threshold = Threshold_ratio*double(max(H(:)));
num_nb = 1;
az_ind=[];el_ind=[];xp_ind=[];yp_ind=[];h_ind=[];
if max(H(:))~=0 && max(H(:))~=1 %scattered segments are noise, skip while loop
    while max(Htemp(:))>= H_Threshold
        [az_ind_temp,el_ind_temp,xp_ind_temp,yp_ind_temp]=ind2sub(size(H),find(Htemp==max(Htemp(:)),1));
        [az_ind_nb,el_ind_nb,xp_ind_nb,yp_ind_nb] = ndgrid(findWrappedInds(az_ind_temp, length(az_bin)-1, num_nb),...
            findWrappedInds(el_ind_temp, length(el_bin)-1, num_nb),...
            max(1,xp_ind_temp-num_nb):min(length(xp_bin)-1,xp_ind_temp+num_nb),...
            findWrappedInds(yp_ind_temp, length(yp_bin)-1, num_nb));
        nbs = sub2ind(size(H),az_ind_nb(:),el_ind_nb(:),xp_ind_nb(:),yp_ind_nb(:));
        Htemp(nbs) = 0;
        az_ind=[az_ind;az_ind_temp];el_ind=[el_ind;el_ind_temp];
        xp_ind=[xp_ind;xp_ind_temp];yp_ind=[yp_ind;yp_ind_temp];
        h_ind=[h_ind;max(Htemp(:))];
        if length(az_ind) > numpeaks 
            break; %select "numpeaks" max H
        end
    end
end

lines = struct('point1',[],'point2',[],'az',[],'el',[],'xp',[],'yp',[],'H',[]);
for l = 1:length(az_ind)
    PtPairInd = H_ptPairInd{az_ind(l),el_ind(l),xp_ind(l),yp_ind(l)}';
    AllPts = indcomb(PtPairInd,:);
    % find the most distant pair
    dists=DistanceBtPoints(currPtCloud(AllPts(:,1),:),currPtCloud(AllPts(:,2),:));
    [len,ii] = max(dists);
    if len < minLength
       continue;
    end       
    lines(l).az = 0.5*(az_bin(az_ind(l)) + az_bin(az_ind(l)+1));
    lines(l).el = 0.5*(el_bin(el_ind(l)) + el_bin(el_ind(l)+1));
    lines(l).xp = 0.5*(xp_bin(xp_ind(l)) + xp_bin(xp_ind(l)+1));
    lines(l).yp = 0.5*(yp_bin(yp_ind(l)) + yp_bin(yp_ind(l)+1));
    lines(l).point1 = currPtCloud(AllPts(ii,1),:);
    lines(l).point2 = currPtCloud(AllPts(ii,2),:);
    lines(l).H = h_ind(l);
end
lines_all = [lines_all; lines(:)];
end

nonemptyIndex = arrayfun(@(array) ~isempty(array.point1),lines_all);
if ~any(nonemptyIndex)
    H=[];az_bin=[];el_bin=[];xp_bin=[];yp_bin=[];Lines = [];
    return;
end
lines_all = lines_all(nonemptyIndex);

if length(lines_all)>3
    [~, ord] = sort([lines_all(:).H],'descend');
    Lines = lines_all(ord(1:3));
else
    Lines = lines_all;
end






