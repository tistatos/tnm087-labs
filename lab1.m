%%
% Lab 1 - HDR
%%
clear

exposureRatio = 2;

load('gfun.mat');

%% plot g function
subplot(2,1,1);
plot(gfun);
subplot(2,1,2);
plot(2.^gfun);


%% read in pictures
for i=1:14
    pictures(:,:,:,i) = imread(strcat('Img',num2str(i),'.tiff'));
    gPictures(:,:,i) = rgb2gray(pictures(:,:,:,i));
end


hiIntVal = max(max(gPictures(:,:,1)));
[hiIntX, hiIntY] = find (gPictures(:,:,1) == hiIntVal);
loIntVal = min(min(gPictures(:,:,14)));
[loIntX, loIntY] = find (gPictures(:,:,14) == loIntVal);
medianGray = median(median(gPictures(:,:,9)));
[medianGrayX, medianGrayY] = find(gPictures(:,:,9) == medianGray);


for i=1:14
    hiIntValues(:,i) = gPictures(hiIntX, hiIntY, i);
    loIntValues(:,i) = gPictures(loIntX(1,1), loIntY(1,1), i);
    medianValues(:,i) = gPictures(medianGrayX(1,1), medianGrayY(1,1), i);
end

plot(hiIntValues, 'r');
hold on
plot(loIntValues, 'g');
plot(medianValues, 'b');