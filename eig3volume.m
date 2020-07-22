function [Lambda1,Lambda2,Lambda3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz)
% matlab code of eig3volume
% calculates the eigen values from the hessian matrix, sorted by abs value. 
% do not output eigenvector
% reference formula in Oliver K Smith:
% Communications of the ACM April 1961 https://doi.org/10.1145/355578.366316
% and "A Method for Fast Diagonalization of a 2x2 or 3x3 Real Symmetric
% Matrix" by M.J. Kronenburg
% 
% [Lambda1,Lambda2,Lambda3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz)
%
%      | Dxx  Dxy  Dxz |
% H =  | Dyx  Dyy  Dyz |
%      | Dzx  Dzy  Dzz |
%
% Dxy == Dyx; Dxz == Dzx; Dyz == Dzy
% det(H) = -Dxz.^2.*Dyy + 2.*Dxy.*Dxz.*Dyz - Dxx.*Dyz.^2 - Dxy.^2.*Dzz + Dxx.*Dyy.*Dzz

% characteristic equation lambda^3 - b*lambda^2 + c*lambda + d = 0
b = Dxx + Dyy + Dzz;
c = Dxx.*Dyy + Dxx.*Dzz + Dyy.*Dzz - Dxy.^2 - Dyz.^2 - Dxz.^2;
d = Dxx.*Dyz.^2 + Dyy.*Dxz.^2 + Dzz.*Dxy.^2 - Dxx.*Dyy.*Dzz - 2.*Dxy.*Dxz.*Dyz;
p = b.^2-3.*c;
q = 2.*b.^3 - 9.*b.*c - 27.*d;

discriminant = acos(q./2./sqrt(p.^3));
Lambda1 = 1/3*(b+2.*sqrt(p).*cos(discriminant./3));
Lambda2 = 1/3*(b+2.*sqrt(p).*cos((discriminant+2*pi)./3));
Lambda3 = 1/3*(b+2.*sqrt(p).*cos((discriminant-2*pi)./3));

% complex values caused by rounding in calculation; use matlab function eig
comp_ind = find(abs(q./2./sqrt(p.^3))>1);
for i=1:length(comp_ind)
    Lambdas = eig([Dxx(i) Dxy(i) Dxz(i);Dxy(i)  Dyy(i)  Dyz(i);Dxz(i)  Dyz(i)  Dzz(i)]);
    Lambda1(i) = Lambdas(1);
    Lambda2(i) = Lambdas(1);
    Lambda3(i) = Lambdas(1);
end

ind_AlargerthanC = find(abs(Lambda1)>abs(Lambda3));
temp = Lambda1(ind_AlargerthanC);
Lambda1(ind_AlargerthanC) = Lambda3(ind_AlargerthanC);
Lambda3(ind_AlargerthanC) = temp;

ind_AlargerthanB = find(abs(Lambda1)>abs(Lambda2));
temp = Lambda1(ind_AlargerthanB);
Lambda1(ind_AlargerthanB) = Lambda2(ind_AlargerthanB);
Lambda2(ind_AlargerthanB) = temp;

ind_BlargerthanC = find(abs(Lambda2)>abs(Lambda3));
temp = Lambda2(ind_BlargerthanC);
Lambda2(ind_BlargerthanC) = Lambda3(ind_BlargerthanC);
Lambda3(ind_BlargerthanC) = temp;






