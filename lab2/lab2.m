%%
% Lab 2 - Operations in the image domain
%%

% Clear any saved data
clear

%% Exercise 1 - part a
N = 512;
M = 0:100;
Canon = imread('CWhite1.jpg');
Holga = imread('HWhite1.jpg');

Canon = imresize(Canon, 0.5);
Holga = imresize(Holga, 0.5);

imshow(Canon);
figure
imshow(Holga);

[X,Y] = meshgrid((1:N));
[T,R] = cart2pol(X-N/2,Y-N/2);

plot(R(N/2,:));
imshow(R,[]);




%% Exercise 1 - part b
SR = R/255;

SR = SR/1.4198;
SR = SR*100;
QR = int16(round(SR));

canonMaskAvg = zeros(100, 1);
holgaMaskAvg = zeros(100, 1);
for m = 1:100
    Maskm = QR == m;
    %Sum over the pixel values
    masksum = sum(sum(Maskm));
    %object point?
    Canonmask = Canon(Maskm == 1);
    Holgamask = Holga(Maskm == 1);
    canonMaskAvg(m) = mean(Canonmask);
    holgaMaskAvg(m) = mean(Holgamask);
    
    
end



plot(canonMaskAvg);
plot(holgaMaskAvg);

canonMaskAvg = canonMaskAvg/max(canonMaskAvg);
%holgaMaskAvg = holgaMaskAvg/max(holgaMaskAvg);

plot(canonMaskAvg);
%plot(holgaMaskAvg);


%% Exercise 2

RedEyes = imread('BoldRedEye.jpg');
load('RedEyeMask.mat');
RedChannel = RedEyes(:,:,1);
imshow(RedChannel);

RedChannel = double(RedChannel)/255;
squareFilter = ones(64);

MFilterImage = imfilter(RedChannel, squareFilter);
EyeFilterImage = imfilter(RedChannel, RedEyeMask);

combinedFilterImage = EyeFilterImage./MFilterImage;

imshow(EyeFilterImage/max(max(EyeFilterImage)));
imshow(MFilterImage/max(max(MFilterImage)));

imshow(combinedFilterImage/max(max(combinedFilterImage)));

quantvalues = quantile(quantile(combinedFilterImage, 0.985), 0.985);

combinedFilterImage = combinedFilterImage >= quantvalues;
BW = imregionalmax(combinedFilterImage);
%imshowpair(BW, RedEyes);
imshow(BW);

%% Exercise 3
GIM = imread('GCPins512.jpg');
HIM = imread('GHPins512.jpg');

imshow(GIM);
[GX, GY] = ginput(3);
for k= 1:3
    text(GX(k),GY(k),num2str(k));
end
hold on
figure
imshow(HIM);
[HX, HY] = ginput(3);

GC = [GX'; GY'; 1 1 1;]';
HC = [HX'; HY'; 1 1 1;]';

x = mldivide(GC,HC);
c = affine2d(x);

Rcb = imref2d(size(GIM))

transFormedImage = imwarp(GIM, c, 'OutputView', Rcb);

imshowpair(transFormedImage, HIM);


