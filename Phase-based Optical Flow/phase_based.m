%% Extract dataset
clear;clc;
% Add path to directory containing 'phase_based.m' and its Utility files
% Utility files - 'phase_based_util.m' and 'refFrame.m'
% Path file to the video
vr = VideoReader('Grooving_5406fps.mp4');
nF = vr.NumberOfFrames; % number of Frames
disp('Finished loading video');
%% Create Gabor filter bank
orientation = [0 90]; % List of orientations for the filter bank
N = length(orientation); % number of filters to be created
wavelength = repmat(10,[1,N]); % set wavelength in pixels per cycle for each filter
bandwidth = repmat(3,[1,N]); % set bandwidth in octaves for each filter
gaborBank = gabor(wavelength,orientation,'SpatialFrequencyBandwidth',bandwidth); % generate filter bank
fprintf('Gabor filter size: (%d,%d)\n',size(gaborBank(:,1).SpatialKernel));
%% Calculate motion
clc; close all;
range1 = 143; % start frame for displacement calculation
range2 = 243; % end frame for displacement calculation
img1 = rgb2gray(vr.read(range1)); % read first frame in the range
horizontal = 1; vertical = 0; % set the orientation of velocity required. Selection both together also valid
sigma = 1; % Standard deviation parameter of gaussian filter used for smoothing image before gabor filter application
do_dilation = 0; % Set it to 1 to dilate the threshold mask. Dilation is a morphological operation and thickens the mask

% Call reference frame utility function. It takes the first frame in the
% range and creates threshold masks and calculates phase information from
% the output of the the gabor filter. All information is stored as object
% variable for faster computation
img1_class = refFrame(img1(370:419,740:789), gaborBank, sigma, do_dilation);

u = []; % stores displacement matrix calculated
tic % for calculating time elapsed between execution of tic and toc
n = range2-range1-1; % number of frames for displacement calculation

f = waitbar(0, 'Starting'); % create a waitbar

for i=range1+1:range2
    img2 = rgb2gray(vr.read(i)); % read frame
    % Call Utility class. Takes reference class object, current image,
    % gabor filters, standard deviation for gaussian smoothing, and
    % orientation of displacements required. Returns displacements in
    % horizontal and vertical directions of each pixel in the image.
    [temp_u, ~] = phase_based_util(img1_class, img2(370:419,740:789), horizontal, vertical,...
        gaborBank, sigma);
    u = cat(3,u, temp_u); % concatenate displacement matrix to vector of displacement matrix
    waitbar((i-range1-1)/n, f, sprintf('Progress: %d %%\n', floor((i-range1-1)/n*100))); % update waitbar
end
toc % end runtime calculation
waitbar(1,f,'Done'); 
close(f); % close waitbar

% imshow(img1_class.thres_mask_h,[]); % Visualize threshold mask
% saveas(gcf,"Phase based image",'fig'); % Save phase based figure

X=[]; % average displacement for each frame
for i = 1:size(u,3)
    X(i) = mean(reshape(u(:,:,i),[],1)); % take average of each displacement matrix frame for scalar displacement
end

fs = 2142; % frame rate of video recording
X = X-mean(X); % make mean zero
disp = X;
X = [X, zeros(1,200)]; % add zeros at the end to calculate smoother fft
f = fs*linspace(0,1,length(X)); % frequency vector
time = (1/fs)*(0:length(disp)-1); % time vector

pixel_res = 83; % microns per pixel 

% Create figure
figure(22);
subplot(211)
plot(time, pixel_res*349*disp/1000,'b', 'LineWidth',2.0); % scale displacement from phase based
xlim([0, size(u,3)/fs]);
ylim([-0.007, 0.007]);
xlabel('Time (s)','Interpreter','latex');
title("Phase based optical flow displacement",...
    'interpreter','latex');
ylabel('Displacement (mm)','Interpreter','latex');
set(gca,'Fontsize',14);
grid minor

subplot(212)
fft_X = (abs(fft(pixel_res*349*disp/1000,length(X))));
plot(f,(fft_X),'b','linewidth',2.0)
ylabel("Magnitude", ...
    'Interpreter','latex')
xlabel("Frequency (Hz)",...
    'Interpreter','latex')
grid minor
xlim([1,fs/2])
set(gca,"FontSize",14,'yscale','linear','xscale','linear')
title('Fast Fourier Transform','Interpreter','latex');
%% Plot region of interest for displacement calcuation
imshow(img1,[])
hold on;
rectangle('Position',[740,370,50,50],...
    'Curvature',[0,0],...
    'LineWidth',4,'LineStyle','-','edgecolor','r')
%% Visualize displacement matrix
close all;
for i = 1:size(u,3)
    imagesc(tanh(5*u(:,:,i)),[-0.18 0.01]);
    pause(0.05);
end
%% Save displacement vector
save('may5_end_mill_2142fps_disp_PHASE_BASED_Box_A.mat','disp');
fprintf('Saved\n');