clear
clc


%% 

[data,dataname]=xlsread('juliendata.xlsx');
imds = imageDatastore('D:\summerschool\figFT\fig0\figFT',...  
    'IncludeSubfolders',true,...  
    'LabelSource','foldernames');  % path to training set

% a=double(imds.Labels);
% labels=de2bi(a)
%% show Labels and Count
data_disp = countEachLabel(imds);
disp(data_disp);
  
%%   get features one spectrum
%   features depend on size and hog
image1 = readimage(imds,1);  %convert to img matrix
features=getfeature(image1);
%% 
  
% get features from the spectrum  
numImages = length(imds.Files); 
timespan=0.5;
cuttimes=60/timespan;
features = zeros(numImages*cuttimes,size(features,2),'single'); % featuresTrain is the single precision
featureslabels=cell(numImages*cuttimes,1);
target=cell(numImages*cuttimes,1);
for i = 1:numImages 
    allimage = readimage(imds,i);  
    for j= 1:cuttimes
        n=325/cuttimes;
        m=(j-1)*n;
        imgcut1=imcrop(allimage,[m 0 n 343]);
        features((i-1)*cuttimes+j,:)=getfeature(allimage);
        target((i-1)*cuttimes+j,:)=dataname(i,6);
        featureslabels((i-1)*cuttimes+j,:)=dataname(i,7);
    end
end  

trainLabels = categorical(featureslabels);  
rng(1);  
t=templateSVM('Standardize',1);
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
title('0.5s crossvaldition accuary (percentage)')
saveas(gcf,'f1per.jpg')
figure(2)
confusionchart(Z,'RowSummary','absolute','Normalization','absolute')
title('0.5s crossvaldition accuary ')
saveas(gcf,'f1.jpg')
kloss=kfoldLoss(cvmdl);

save('cvmdl');
save('classifier');
save('target');
save('features');
save('featureslabels');
save('kloss');



%%

function feature=getfeature(RGB)
I=rgb2gray(RGB);
BW = imbinarize(I);
B=bwboundaries(BW,'noholes');
maxf=zeros(length(B),1);
minf=zeros(length(B),1);
for k = 1:length(B)
   boundary = B{k};
   maxf(k)=max(boundary(:,1));
   minf(k)=min(boundary(:,1));
end
maxf=max(maxf);
minf=min(minf);
fspan=maxf-minf;
m=mean(I,'all');
pd=bwarea(BW);
kurtosisBW=kurtosis(RGB,1,'all');
Skewness = skewness(RGB,1,'all');
feature=[maxf,minf,fspan,m,pd,kurtosisBW,Skewness];
end

