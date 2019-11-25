function principalCurve=computePC(X,l_KR,gridPtUnit)
% approximates principal curve from points
% gridPtUnit is in the unit of pixels

Xn = X - mean(X); % move data to the center
coeff = pca(Xn); 
Xtr=(coeff\(Xn'))'; % transfer to PCA-x, PCA-y coordinates
gridPt=min(Xtr(:,1)):gridPtUnit:max(Xtr(:,1));
Xtr_sort = sortrows(Xtr); % sort projection points

D=repmat(gridPt,size(X,1),1)-repmat(Xtr_sort(:,1),1,size(gridPt,2));
K=exp(-D.^2/(2*l_KR^2)); % kernel
Kn=K./repmat(sum(K),size(X,1),1); % normalize
Str=[gridPt' Kn'*Xtr_sort(:,2)];
principalCurve = (coeff*Str')'+mean(X);
