RGB=imread('127Moving arm faster towards radar, slower awayJulien.png');
I=rgb2gray(RGB);
imshow(I);
BW = imbinarize(I);
Peak2Peak = peak2peak(RGB)