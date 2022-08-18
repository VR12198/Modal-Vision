function edge_final = edge_sobel(img,Th)
img = double (img);
%Value for Thresholding

%Filter for horizontal and vertical direction
KGx = [1, 0, -1; 2, 0, -2; 1, 0, -1];
KGy = [1, 2, 1; 0, 0, 0; -1, -2, -1];

%Convolution by image by horizontal and vertical filter
Filtered_X = conv2(img, KGx, 'same');
Filtered_Y = conv2(img, KGy, 'same');
pan=size(img,1);
leb=size(img,2);
%Calculate magnitude
magnitude = (Filtered_X.^2) + (Filtered_Y.^2);
BW = sqrt(magnitude);

%Hysteresis Thresholding
Th = Th * max(max(BW));
T_res = zeros (pan, leb);
for i = 1  : pan
    for j = 1 : leb
        if (BW(i, j) < Th)
            T_res(i, j) = 0;
        else
            T_res(i, j) = 1;
        %Using 8-connected components
        end
    end
end
T_res(1:5,:) = 0;
T_res(:,1:5) = 0;
T_res(end-5:end,:) = 0;
T_res(:,end-5:end) = 0;
edge_final = T_res;

%% Show final edge detection result
% close all
% figure
% imshow(imdilate(edge_final,strel('disk',2)));