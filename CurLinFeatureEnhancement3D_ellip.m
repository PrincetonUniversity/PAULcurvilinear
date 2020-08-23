function outputIM=CurLinFeatureEnhancement3D_ellip(contrastEnhancedIM, PSFsigma, pixelSize, Alpha, Beta, cRatio)

% Modified from "FrangiFilter3D" by Marc Schrijver and Dirk-Jan Kroon (University of Twente)
% https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter

scale2PSFratio = [2 2 2];
TwoAlphaSq = 2*Alpha^2;
TwoBetaSq  = 2*Beta^2;

% format input
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

%%%%%% fine screening for different pixelsize, PSFsigma in x, y, z %%%%%%%
% [sigmasX,sigmasY,sigmasZ] = meshgrid(...
%     floor(scale2PSFratio*PSFsigma(1)/pixelSize(1)):...
%     0.25:ceil(scale2PSFratio*PSFsigma(1)/pixelSize(1)),...
%     floor(scale2PSFratio*PSFsigma(2)/pixelSize(2)):...
%     0.25:ceil(scale2PSFratio*PSFsigma(2)/pixelSize(2)),...
%     floor(scale2PSFratio*PSFsigma(3)/pixelSize(3)):...
%     0.25*round(scale2PSFratio*PSFsigma(3)/pixelSize(3)):ceil(scale2PSFratio*PSFsigma(3)/pixelSize(3)));
% 
% sigmasX = sigmasX(:);
% sigmasY = sigmasY(:);
% sigmasZ = sigmasZ(:);
%%%%%% fine screening for different pixelsize, PSFsigma in x, y, z %%%%%%%

%%%%%% quick "symmetric" screening to save time %%%%%%%%
[sigmasXY,sigmasZ] = meshgrid(...
    floor(scale2PSFratio(1)*PSFsigma(1)/pixelSize(1)):...
    0.25:ceil(scale2PSFratio(1)*PSFsigma(1)/pixelSize(1)),...
    floor(scale2PSFratio(3)*PSFsigma(3)/pixelSize(3)):...
    0.25:ceil(scale2PSFratio(3)*PSFsigma(3)/pixelSize(3)));
sigmasX = sigmasXY(:);
sigmasY = sigmasXY(:);
sigmasZ = sigmasZ(:);
%%%%%% quick "symmetric" screening to save time %%%%%%%%

%%%%%% quick "asymmetric" screening to save time %%%%%%%%
% [sigmasX,sigmasY] = meshgrid(...
%     floor(scale2PSFratio*PSFsigma(1)/pixelSize(1)):...
%     0.25:ceil(scale2PSFratio*PSFsigma(1)/pixelSize(1)),...
%     floor(scale2PSFratio*PSFsigma(2)/pixelSize(2)):...
%     0.25:ceil(scale2PSFratio*PSFsigma(2)/pixelSize(2)));
% sigmasZ = repmat(scale2PSFratio*PSFsigma(3)/pixelSize(3), size(sigmasX));
% sigmasX = sigmasX(:);
% sigmasY = sigmasY(:);
% sigmasZ = sigmasZ(:); 
%%%%%% quick "asymmetric" screening to save time %%%%%%%%

% Frangi filter for all sigmas
for i = 1:length(sigmasX)
    
    %     z2xy_scaleRatio = 0.23/0.076;
    % Calculate 3D hessian
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(contrastEnhancedIM,[sigmasX(i),sigmasY(i),sigmasZ(i)]);
    
    % Correct for scaling
    Dxx = (sigmasX(i)^2)*Dxx;
    Dyy = (sigmasY(i)^2)*Dyy;
    Dzz = (sigmasZ(i)^2)*Dzz;
    Dxy = (sigmasX(i)*sigmasY(i))*Dxy;
    Dxz = (sigmasX(i)*sigmasZ(i))*Dxz;
    Dyz = (sigmasY(i)*sigmasZ(i))*Dyz;
    
    % Calculate eigen values
    [Lambda1,Lambda2,Lambda3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    clear Dxx Dyy Dzz Dxy Dxz Dyz;
    
    % Calculate absolute values of eigen values
    LambdaAbs1=abs(Lambda1);
    LambdaAbs2=abs(Lambda2);
    LambdaAbs3=abs(Lambda3);
    
    % The Vesselness Features
    Ra=LambdaAbs2./LambdaAbs3;
    Rb=LambdaAbs1./sqrt(LambdaAbs2.*LambdaAbs3);
    
    % Second order structureness. S = sqrt(sum(L^2[i])) met i =< D
    S = sqrt(LambdaAbs1.^2+LambdaAbs2.^2+LambdaAbs3.^2);
    TwoCSq=2.*(cRatio*max(S(:))).^2;% cRatio's maximum Hessian F-norm
    clear LambdaAbs1 LambdaAbs2 LambdaAbs3
    
    %Compute Vesselness function
    expRa = 1-exp(-Ra.^2./TwoAlphaSq);
    expRb =    exp(-Rb.^2./TwoBetaSq);
    expS  = 1-exp(-S.^2./TwoCSq);
    clear S A B C Ra Rb
    
    %Compute Vesselness function
    Voxel_data = expRa.* expRb.* expS;
    clear expRa expRb expRc;
    Voxel_data(Lambda2 > 0)=0; Voxel_data(Lambda3 > 0)=0;
    
    % Remove NaN values
    Voxel_data(~isfinite(Voxel_data))=0;
    Voxel_data(isnan(Voxel_data))=0;
    
    % keep the largest volume when screening sigmas
    if i==1
        outputIM = Voxel_data;
    else
        outputIM=max(outputIM,Voxel_data);
    end
end

