function [f,H] = addNoise(Img)
%Prev_f = im2double(rgb2gray(imread(path)));
%Generate the blurring filter
hsize = 5;
std = 7; % 5 %
H = fspecial('gaussian',hsize,std);
Bim = imfilter(Img,H,'circular','conv');
sigma = 10/255;
f = imnoise(Bim,'gaussian',0,sigma*sigma);
%f = imnoise(Bim,'salt & pepper');
%pathT = erase(path,".bmp");
%NewPath = strcat(pathT,"Noise.jpg");
%imwrite(f,NewPath);
end