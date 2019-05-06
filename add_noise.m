clc;
clear;
close all;

image=imread('rect_full.jpg');
im=rgb2gray(image);    
figure;
imshow(mat2gray(im));
title('reference image');

%%%%%%%%%%%% warping %%%%%%%%%%%%%%
warp=[1 0.5 -20; 0 1 0; 0 0 1];
[A,B]=size(im);
nx=1:B;
ny=1:A;
wim = iat_inverse_warping(im, warp, 'affine', nx, ny); 
figure;
imshow(mat2gray(wim));


% %%%%%%%%%%% noise %%%%%%%%%%%%%
% image=imnoise(image,'gaussian',0,0.02);
% % mask=fspecial('gaussian',[5,5],1);
% % image=imfilter(image,mask);
% figure;
% imshow(image);

%%%%%%%%% imwrite %%%%%%%%%%%%%
imwrite(mat2gray(wim),'rect_sheer.jpg');

% %%%%%%%%%% corner %%%%%%%%%%%%%
% C1=detectHarrisFeatures(image);
% C1=C1.selectStrongest(4);
% hold on;
% plot(C1);
% 
% %%%%%%%%%% warped image %%%%%%%%%%%%
% w=imread('04.png');
% w=rgb2gray(w);
% figure;
% imshow(w);
% title('warped image');
% C2=detectHarrisFeatures(w);
% C2=C2.selectStrongest(9);
% hold on;
% plot(C2);
% 
% %%%%%%%%%% matching features %%%%%%%%%%%
% %%Extract the neighborhood features.
% [features1,valid_points1] = extractFeatures(image,C1);
% [features2,valid_points2] = extractFeatures(w,C2);
% 
% %%Match the features.
% indexPairs = matchFeatures(features1,features2);
% 
% %%Retrieve the locations of the corresponding points for each image.
% matchedPoints1 = valid_points1(indexPairs(:,1),:);
% matchedPoints2 = valid_points2(indexPairs(:,2),:);
% 
% %%Visualize the corresponding points. You can see the effect of translation 
% %%between the two images despite several erroneous matches.
% figure; 
% showMatchedFeatures(image,w,matchedPoints1,matchedPoints2);
% title('matching image');
% legend('matched points 1','matched points 2');
% 
% % features1.Features(1,:)


