%% Extract dataset
% Code motivated from work by Pulkit Gupta (https://doi.org/10.1016/j.cirpj.2020.11.012)
clc; clear; close;
% Add path to directory containing 'optical_flow.m' and its Utility files
% Utility files - 'HS2.m', 'edge_sobel.m', 'edge_canny.m' and 'computeDerivates.m'
% Path file to the video
vr = VideoReader('C:\Users\Varun\Desktop\Lab\DataSet\stiff_endmill_lens24_fps_2142_clipped_video.mp4');
nF = vr.NumberOfFrames; % number of Frames
X = zeros(nF,1); % displacement vector
disp('Finished loading video');
%%
clc; close all;
range1 = 144; % start frame for displacement calculation
range2 = 243; % end frame for displacement calculation
sigma = 1; % Standard deviation parameter of gaussian filter used for smoothing image before gabor filter application
methods = 'sobel';% {'HS';'sobel';'prewitt'};
edge_type = 1; % {0 or 1};

img1 = rgb2gray(vr.read(range1));
img1 = img1(370:419,740:789); % Crop first frame to region of interest
n = range2-range1-1; % Number of frames for displacement calculation
f = waitbar(0, 'Starting'); % initialize waitbar
for i=range1+1:range2
    img2 = rgb2gray(vr.read(i));
    img2 = img2(370:419,740:789); % Crop frames to region of interest
    s = HS2(img1,img2,sigma,methods); % caluculate full field displacement
    BW = edge_sobel(img2,0.37); % Detect edge to get edge mask, adjust threshold to find tune sobel mask
%     BW = imdilate(BW, strel('disk',1)); % perform mask dilation. Dilation is a morphological operation and thickens the mask
    if(edge_type)
        edge_req = s.*BW; % Mask the displacement
    else
        edge_req = s;
    end
    if nnz(edge_req)~=0
        X(i) = sum(edge_req(:))/nnz(edge_req); % Calculate mean of displacement values within the mask
    end
    waitbar((i-range1-1)/n, f, sprintf('Progress: %d %%\n', floor((i-range1-1)/n*100))); % update waitbar
end
waitbar(1,f,'Done');
close(f) % close waitbar
X = X(range1:range2);

fs = 2142; % camera fps
X = (X - mean(X)); % make mean zero
f = fs*linspace(0,1,length(X)); % frequency vector
time = (1/fs)*(0:length(X)-1); % time vector
pixel_res = 83; % microns per pixel  
disp = X*pixel_res/1000; % scaled displacement in mm

%Plot figures
figure()
subplot(2,1,1)
plot(time,disp,'-r','Linewidth',2)
xlim([min(time),max(time)]);
ylim([-0.007, 0.007]);
xlabel('time [s]','Interpreter','latex');
ylabel('Displacement (mm)','Interpreter','latex')
grid minor;
title('Displacement','Interpreter','latex')
set(gca,'Fontsize',14);
subplot(2,1,2)
plot(f,(abs(fft(disp))),'r','linewidth',2.0)
xlim([0 fs/2])
ylabel("Magnitude", ...
    'Interpreter','latex')
xlabel("Frequency (Hz)",...
    'Interpreter','latex')
grid minor
set(gca,"FontSize",14)
title('Fast Fourier Transform','Interpreter','latex');
hold on

clc;
%% SAVE DISPLACEMENT VECTOR
save('may5_end_mill_2142fps_disp_INTENSITY_BASED_Box_A.mat','disp');
fprintf('Saved\n');