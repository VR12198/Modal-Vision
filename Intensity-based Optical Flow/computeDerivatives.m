function [fx, fy] = computeDerivatives(im1, im2, methods)

if size(im2,1)==0
    im2=zeros(size(im1));
end

if(methods==1)
    % Horn-Schunck original method
    fx = conv2(im1,0.25* [-1 1; -1 1],'same') + conv2(im2, 0.25*[-1 1; -1 1],'same');
    fy = conv2(im1, 0.25*[-1 -1; 1 1], 'same') + conv2(im2, 0.25*[-1 -1; 1 1], 'same');
elseif(methods==2)
    % Sobel
    fx = conv2(im1,(1/16)*[1 2 1; 0 0 0; -1 -2 -1],'same') + conv2(im2,(1/16)*[1 2 1; 0 0 0; -1 -2 -1],'same');
    fy = conv2(im1,(1/16)*[-1 0 1; -2 0 2; -1 0 1],'same') + conv2(im2,(1/16)*[-1 0 1; -2 0 2; -1 0 1],'same');
else
    % Prewitt
    fx = conv2(im1,(1/12)*[-1 -1 -1; 0 0 0; 1 1 1],'same') + conv2(im2,(1/12)*[-1 -1 -1; 0 0 0; 1 1 1],'same');
    fy = conv2(im1,(1/12)*[-1 0 1; -1 0 1; -1 0 1],'same') + conv2(im2,(1/12)*[-1 0 1; -1 0 1; -1 0 1],'same');
end