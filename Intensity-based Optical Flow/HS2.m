function [s] = HS2(im1, im2, sigma, methods)

if size(size(im1),2)==3
    im1=rgb2gray(im1);
end
if size(size(im2),2)==3
    im2=rgb2gray(im2);
end
im1=double(im1);
im2=double(im2);

hsize = 6*sigma;
if(mod(hsize,2)==0)
    hsize = hsize+1;
end

if(sigma~=0)
    B = fspecial('gaussian',hsize,sigma);
    im1 = conv2(im1, B, 'same');
    im2 = conv2(im2, B, 'same');
end

%% Estimate spatiotemporal derivatives
[fx, fy] = computeDerivatives(im1, im2, methods);
s = (im2-im1)./((fx.^2+fy.^2).^0.5);
s(isinf(s)|isnan(s))=0;