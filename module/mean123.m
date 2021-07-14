RGB=imread('1WalkingAleksandar.png');
I=rgb2gray(RGB);
imshow(I);
BW = imbinarize(I);
m=mean(I,'all')