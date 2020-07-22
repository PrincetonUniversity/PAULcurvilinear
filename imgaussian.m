function I=imgaussian(I,sigma,siz)
% modified to allow sigma to be 1x2 or 1x3 vectors
% so that sigma is not equal at different axes
% x: 1st dimension; y: 2nd dimension; z: 3rd dimension

% IMGAUSSIAN filters an 1D, 2D color/greyscale or 3D image with an
% Gaussian filter. This function uses for filtering IMFILTER or if
% compiled the fast  mex code imgaussian.c . Instead of using a
% multidimensional gaussian kernel, it uses the fact that a Gaussian
% filter can be separated in 1D gaussian kernels.
%
% J=IMGAUSSIAN(I,SIGMA,SIZE)
%
% inputs,
%   I: The 1D, 2D greyscale/color, or 3D input image with
%           data type Single or Double
%   SIGMA: The sigma used for the Gaussian kernel
%   SIZE: Kernel size (single value) (default: sigma*6)
%
% outputs,
%   J: The gaussian filtered image
%
% note, compile the code with: mex imgaussian.c -v
%
% example,
%   I = im2double(imread('peppers.png'));
%   figure, imshow(imgaussian(I,10));
%
% Function is written by D.Kroon University of Twente (September 2009)

if isscalar(sigma)
    sigma = [sigma sigma sigma];
elseif length(sigma) == 2
    sigma = [sigma(1) sigma(1) sigma(2)];
end

if(~exist('siz','var')), siz=sigma.*6; end

% Make 1D Gaussian kernels
x=-ceil(siz(1)/2):ceil(siz(1)/2);Hx = exp(-(x.^2/(2*sigma(1)^2)));
y=-ceil(siz(2)/2):ceil(siz(2)/2);Hy = exp(-(y.^2/(2*sigma(2)^2)));
z=-ceil(siz(3)/2):ceil(siz(3)/2);Hz = exp(-(z.^2/(2*sigma(3)^2)));
Hx = Hx/sum(Hx(:));Hy = Hy/sum(Hy(:));Hz = Hz/sum(Hz(:));

% Filter each dimension with the 1D Gaussian kernels
Hx=reshape(Hx,[length(Hx) 1 1]);
Hy=reshape(Hy,[1 length(Hy) 1]);
Hz=reshape(Hz,[1 1 length(Hz)]);
I=imfilter(imfilter(imfilter(I,Hx, 'same','replicate'),Hy, 'same','replicate'),Hz, 'same','replicate');

