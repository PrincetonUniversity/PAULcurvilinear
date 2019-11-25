function y_max=findMaxthruInterpolation (x,y)
% modified from M Arthington's splineMaximaMinima (09/12/2009)
% https://www.mathworks.com/matlabcentral/fileexchange/26144-local-maxima-and-minima-of-a-pp-spline
% make discriminant strictly > 0 
% compare all local maxima with all grid values, select the largest

ppSpline = spline(x(:),y(:));
[breaks,coefs]=unmkpp(ppSpline);
discriminant = (4*coefs(:,2).^2-12.*coefs(:,1).*coefs(:,3));
logicalPiecesWithStationary = discriminant>0; % "=0": single inflection point; "<0" monotomic -> changed to ">" only
t2 = breaks(logicalPiecesWithStationary)' + (-2*coefs(logicalPiecesWithStationary,2)-sqrt(discriminant(logicalPiecesWithStationary)))./(6*coefs(logicalPiecesWithStationary,1));
t2IsStationary = t2<breaks([false; logicalPiecesWithStationary])' & t2>breaks(logicalPiecesWithStationary)';
y_max = max(ppval(ppSpline,[breaks(1);breaks(end);t2(t2IsStationary)]));	


