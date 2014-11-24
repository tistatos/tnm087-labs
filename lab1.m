%%
% Lab 1 - HDR
%%

% Clear any saved data
clear

% exposure ratio between images
exposureRatio = 2;
% load g function = log2(f^-1(z_ij))
load('gfun.mat');

% plot g function
subplot(2,1,1);
plot(gfun);
subplot(2,1,2);
plot(2.^gfun);


% read picture into an array
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

% Plot hi,low and median
figure
plot(hiIntValues, 'r')
hold on
plot(medianValues, 'g')
plot(loIntValues, 'b')
%%
figure
tic

%calculate weight image (could be done with one pic-loop but it doesnt seem to work)
weight = pictures;
for pic=1:14
    for y=1:683
        for x=1:1024
            for c=1:3

                value = pictures(y,x,c,pic);
                %weight if value is greater than 128
                if(value > 128)
                    value = 255-value;
                end
                weight(y,x,c,pic) = value;
            end
        end
    end
end

%summarize all pictures in to one weight picture
weightPic = weight(:,:,:,1);
for pic=2:14
   weightPic = weightPic + weight(:,:,:,pic);
end

%normalize weight values
weightPic= double(weightPic)/(255);

%Get irradiance values from gfun and weight it
irradiancePictures = double(pictures);
finv= (2.^gfun);
for pic=1:14
    value = pictures(:,:,:,pic);
    irValue = finv(value+1)/(exposureRatio*pic);
    irradiancePictures(:,:,:,pic) = irValue.*weightPic(:,:,:);
end

%Summarize to finaly picture
finalpic = irradiancePictures(:,:,:,1);
for pic=2:14
   finalpic = imadd(finalpic,irradiancePictures(:,:,:,pic));
end

%tonemap and show the result
imshow(tonemap(finalpic));
toc