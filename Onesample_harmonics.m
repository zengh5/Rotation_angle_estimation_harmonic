
% Clear all variables and close all figures
clear
close all;
% 1. parameter
theta = 5;    % set theta = 5 to reproduce Fig. 3 of the paper
block_size = 128;
% The threshold to determine how much the estimated angle can deviate from
% the ground truth.
threshold = 0.5;
% nearest interpolation introduce stronger harmonics
INT_METHOD = 'nearest';   
% INT_METHOD = 'bilinear';
% INT_METHOD = 'bicubic';

% 2. Image_folder
% Folder = 'H:\Imgdatabase\BossBase-1.0-raw-gray';
Image = 'barbara.bmp';
% 2.1 Read the image
img = imread(Image);
% Important: to avoid the CFA artifacts, for images that directly from the 
% cameras, e.g., images from BossRaw, please uncomment the following line
% img = img(1:2:end,1:2:end);
% 2.2 Rotating
img_rot = imrotate(img, theta, INT_METHOD, 'crop');
% [M, N] = size(img_rot);
% img_rot = img_rot(round(M/2)-41-block_size:round(M/2)+42+block_size, ...
%     round(N/2)-41-block_size:round(N/2)+42+block_size);
figure, imshow(img_rot)

% 3 Angle estimation
% 3.1 Calculating the cyclos spectrum
Txx = CalculateTxx(img_rot, block_size*2);
% 3.2 Remove components near DC
Txx_shift = fftshift(Txx);
[M, N] = size(Txx_shift);
Txx_shift(block_size-1 : block_size+2, block_size-1 : block_size+2) = 0;
Txx = fftshift(Txx_shift);
figure, imshow(Txx,[ ])
% 3.3 Get the estimated rotation angle with the methods of Padin et al. and 
% Chen et al.
estimatedAngle = EstimateAngle(Txx,0);
angle_Chen = estimatedAngle(1);
result_Chen= (abs(angle_Chen-theta)<=threshold);
angle_Padin = estimatedAngle(4);
result_Padin= (abs(angle_Padin-theta)<=threshold);
% 3.4 Get the estimated angle with the proposed method.
% N_har is the highest order of the harmonics to be aggregated
N_har = 2;
estimatedAngle_har = EstimateAngle_har(Txx,0, N_har);
angle_har = estimatedAngle_har;
result_har= (abs(angle_har-theta)<=threshold);

% 4. Show the peaks found
% 4.1 Enhance the spectrum for visualization
Txx_interested = Txx;
Txx_interested(Txx<100)=0;
A = fspecial('average', [7,7]);
Txx_med7_100 = imfilter(Txx_interested, A);
Txx_med = medfilt2(Txx);
Txx_med(Txx_med7_100<=0)=0;
figure, imshow(Txx_med, [], 'border', 'tight')
hold on,
% 4.2 Show the four arcs where Chen's method searches peaks from 
fx = 1: N;
fy1 = ((1-(fx/N).^2).^0.5)*N;
fy2 = N+1-((1-((N+1-fx)/N).^2).^0.5)*N;
fy3 = N+1-((1-(fx/N).^2).^0.5)*N;
fy4 = ((1-((N+1-fx)/N).^2).^0.5)*N;
plot(fx, fy1, 'b--', fx, fy2, 'w--',fx, fy3, 'b--', fx, fy4, 'w--', ...
    'linewidth', 2, 'color', [0.5, 0.5, 0.5])

Done = 1;

