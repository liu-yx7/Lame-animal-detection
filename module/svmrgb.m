clear
clc
[data,dataname]=xlsread('juliendata.xlsx');
imds = imageDatastore('D:\RADAR\RGB\figureRGB_classified',...  
    'IncludeSubfolders',true,...  
    'LabelSource','foldernames');  % path to training set


%% show Labels and Count
data_disp = countEachLabel(imds);
disp(data_disp);
  
%%   get features one spectrum(hog???)
%   features depend on size and hog
imageSize = [256,256]; 
image1 = readimage(imds,1);  %convert to img matrix
image1 = imresize(image1,imageSize);
[features, visualization] = extractHOGFeatures(image1);  %get hog features
  
% get features from the spectrum  
numImages = length(imds.Files);  
features = zeros(numImages,size(features,2),'single'); % featuresTrain is the single precision
featureslabels=cell(numImages,1);
target=cell(numImages,1);
for i = 1:numImages  
    allimage = readimage(imds,i);  
    allimage = imresize(allimage,imageSize);
    features(i,:) = extractHOGFeatures(allimage);
    target(i,:)=dataname(i,6);
    featureslabels(i,:)=dataname(i,7);
end  

trainLabels = categorical(featureslabels);  
rng(1);  
t=templateSVM('Standardize',1)
% start svm training 
classifier = fitcecoc(features,trainLabels,'Learners',t,...
    'ClassNames',{'Bending to pick up and back up','Circling arm forwards','Clapping','Moving arm faster towards radar, slower away','Moving arm slower towards radar, faster away','Sitting and standing','Walking'});  
CVclassifier=crossval(classifier);

cvmdl=CVclassifier;
labels = kfoldPredict(cvmdl);
Y=categorical(trainLabels);
X=categorical(labels);
ConfMat = confusionchart(Y,X);
Z=confusionmat(Y,X);
figure(1);
confusionchart(Z,'RowSummary','row-normalized','Normalization','row-normalized')
title('60s crossvaldition accuary (percentage)')
%saveas(gcf,'f1per.jpg')
figure(2)
confusionchart(Z,'RowSummary','absolute','Normalization','absolute')
title('60s crossvaldition accuary ')
%saveas(gcf,'f1.jpg')

save('cvmdl');