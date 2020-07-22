function [az, el, xp, yp] = convertTo4ParmLineRepDeg_cart(Pt1,Pt2)
% converts a 3D straight line L={r|r=p+t*unitVec{b}} to Roberts' four-parameter
% representation: az [-pi, pi], el [-pi/2, pi/2], rho>0, alpha [-pi, pi]
% Pt1 and Pt2 are two 3D points that define this line
Pt1 = Pt1(:)';
Pt2 = Pt2(:)';

p = Pt1;
if Pt1(3)>=Pt2(3) 
    b = Pt1-Pt2; 
else
    b = Pt2-Pt1; 
end % to ensure bz >= 0
b = b./norm(b);
[az,el,~] = cart2sph(b(1),b(2),b(3));
bx = b(1);
by = b(2);
bz = b(3);

RotMat = [1-bx.^2./(1+bz), -bx.*by./(1+bz), bx;...
    -bx.*by./(1+bz), 1-by.^2./(1+bz),by;...
    -bx, -by, bz]; % RotMat*[0;0;1] = b

xp = p * (RotMat*[1;0;0]); % = rho*cos(alpha)
yp = p * (RotMat*[0;1;0]); % = rho*sin(alpha)
% [alpha,rho] = cart2pol(xp,yp);

az = rad2deg(az);
el = rad2deg(el);
% alpha = rad2deg(alpha);
