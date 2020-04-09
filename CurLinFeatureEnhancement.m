function outputIM=CurLinFeatureEnhancement(contrastEnhancedIM, PSFsigma, pixelSize, Beta, cRatio)

% Modified from "FrangiFilter2D" by Marc Schrijver and Dirk-Jan Kroon (University of Twente)
% https://www.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter

sigmas = floor(2*PSFsigma/pixelSize):0.1:ceil(2*PSFsigma/pixelSize);
TwoBetaSq  = 2*Beta^2;

% Make matrices to store all filterd images
ALLfiltered=zeros([size(contrastEnhancedIM) length(sigmas)]);

% Frangi filter for all sigmas
for i = 1:length(sigmas)
    % Make 2D hessian
    [Dxx,Dxy,Dyy] = Hessian2D(contrastEnhancedIM,sigmas(i));
     
    % Correct for scale
    Dxx = (sigmas(i)^2)*Dxx;
    Dxy = (sigmas(i)^2)*Dxy;
    Dyy = (sigmas(i)^2)*Dyy;
   
    % Calculate (abs sorted) eigenvalues and vectors
    [Lambda1,Lambda2,~,~,~,~]=eig2image(Dxx,Dxy,Dyy);
    
    Rb_sq = (Lambda1./Lambda2).^2;
    S_sq = Lambda1.^2 + Lambda2.^2;
    TwoCSq=cRatio*max(S_sq(:));% cRatio's maximum Hessian F-norm
   
    % Compute the output image
    Ifiltered = exp(-Rb_sq.^2/TwoBetaSq) .*(ones(size(contrastEnhancedIM))-exp(-S_sq/TwoCSq));
    Ifiltered(Lambda2>0)=0;% remove dark tubular background
    % store the results in 3D matrices
    ALLfiltered(:,:,i) = Ifiltered;
end

ALLfiltered(isnan(ALLfiltered)) = 0;

% Return for every pixel the value of the scale(sigma) with the maximum 
% output pixel value

%%%%%%%%%%%% Return maximum on the gridded sigma values (saves time)
% if length(sigmas) > 1
%     [outIm,~] = max(ALLfiltered,[],3);
%     outputIM = reshape(outIm,size(contrastEnhancedIM));
% else
%     outputIM = reshape(ALLfiltered,size(contrastEnhancedIM));
% end

%%%%%%%%%%%% Return the interpolated maximum (not necessarily on the sigma grid)
ALLfiltered_flatten = reshape(ALLfiltered,size(ALLfiltered,1)*...
    size(ALLfiltered,2),size(ALLfiltered,3));
Y = zeros(size(ALLfiltered_flatten,1),1);
parfor y = 1:length(Y)
    Y(y) = findMaxthruInterpolation(sigmas,ALLfiltered_flatten(y,:));
end
% G = (1:size(ALLfiltered_flatten,1))';
% Y0 = splitapply(@(y) findMaxthruInterpolation(sigmas,y),...
%     ALLfiltered_flatten, G);
outputIM = reshape(Y,size(ALLfiltered,1),size(ALLfiltered,2));

