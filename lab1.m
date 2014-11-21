%%
% Lab 1 - HDR
%%
clear

exposureRatio = 2;

load('gfun.mat');

%%% plot g function
subplot(2,1,1);
plot(gfun);
subplot(2,1,2);
plot(2.^gfun);
%%
% read in pictures
for i=1:14
    pictures(:,:,:,i) = imread(strcat('Img',num2str(i),'.tiff'));
    gPictures(:,:,i) = rgb2gray(pictures(:,:,:,i));
end

% get max/median/min values
hiIntVal = max(max(gPictures(:,:,1)));
[hiIntY, hiIntX] = find (gPictures(:,:,1) == hiIntVal);
loIntVal = min(min(gPictures(:,:,14)));
[loIntY, loIntX] = find (gPictures(:,:,14) == loIntVal);
medianGray = median(median(gPictures(:,:,9)));
[medianGrayY, medianGrayX] = find(gPictures(:,:,9) == medianGray);

% get ALL the values!
for i=1:14
    hiIntValues(i) = gPictures(hiIntY, hiIntX, i);
    loIntValues(i) = gPictures(loIntY(1,1), loIntX(1,1), i);
    medianValues(i) = gPictures(medianGrayY(1,1), medianGrayX(1,1), i);
end

%%% Plot hi,low and median
% figure
% plot(hiIntValues, 'r')
% hold on
% plot(medianValues, 'g')
% plot(loIntValues, 'b')
tic


weight = pictures;
for pic=1:14
    for y=1:683
        for x=1:1024
            for c=1:3
                value = pictures(y,x,c,pic);
                if(value > (255/2))
                    value = 255-value;
                end
                weight(y,x,c,pic) = value;
            end
        end
    end
end

weightPic = weight(:,:,:,1);
for pic=2:14
   weightPic = weightPic + weight(:,:,:,pic);
end

weightPic= double(weightPic)/255;

irradiancePictures = double(pictures);
finv= (2.^gfun);
for pic=1:14
    value = pictures(:,:,:,pic);
    irValue = finv(value+1)/(2*(pic));
    irradiancePictures(:,:,:,pic) = irValue.*weightPic(:,:,:);
end

finalpic = irradiancePictures(:,:,:,1);

for pic=2:14
   finalpic = imadd(finalpic,irradiancePictures(:,:,:,pic));
end
imshow(tonemap(finalpic));
toc