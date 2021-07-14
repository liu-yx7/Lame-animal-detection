RGB=imread('2WalkingAleksandar.png');
I=rgb2gray(RGB);
imshow(I);
BW = imbinarize(I);
pd=bwarea(BW)