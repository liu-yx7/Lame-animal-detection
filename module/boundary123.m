RGB=imread('31ClappingFrancesco.png');
I=rgb2gray(RGB);
imshow(I)
BW = imbinarize(I);
[B,L,n,A] = bwboundaries(BW,'noholes');
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