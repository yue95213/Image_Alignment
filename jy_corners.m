 clc;
clear;
close all;

image=imread('03.png');
image=rgb2gray(image);
figure;
imshow(image);
title('reference image');

%%%%%%%%%% corner %%%%%%%%%%%%%
C1=detectHarrisFeatures(image);
hold on;
plot(C1);

%%%%%%%%%% warped image %%%%%%%%%%%%
w=imread('04.png');
w=rgb2gray(w);
figure;
imshow(w);
title('warped image');
C2=detectHarrisFeatures(w);
hold on;
plot(C2);

%%%%%%%%%% matching features %%%%%%%%%%%
%%Extract the neighborhood features.
[features1,valid_points1] = extractFeatures(image,C1);
[features2,valid_points2] = extractFeatures(w,C2);

%%Match the features.
indexPairs = matchFeatures(features1,features2);

%%Retrieve the locations of the corresponding points for each image.
matchedPoints1 = valid_points1(indexPairs(:,1),:);
matchedPoints2 = valid_points2(indexPairs(:,2),:);

%%Visualize the corresponding points. You can see the effect of translation 
%%between the two images despite several erroneous matches.
% figure; 
% showMatchedFeatures(image,w,matchedPoints1,matchedPoints2);
% title('matching image');
% legend('matched points 1','matched points 2');

%%Affine matrix
% transform=[eye(2) zeros(2,1)];
length=length(matchedPoints1);
P1=[matchedPoints1.Location  zeros(length,1)];
P2=[matchedPoints2.Location  zeros(length,1)];


transform=estimateGeometricTransform(matchedPoints1,matchedPoints2,'affine');
figure;
plot(P1(:,1:2),'ro');
hold on;
tform=transform.T(1:2,:);
tform=[tform; 0 0 1];
temp=tform*P1';
temp=temp';
plot(temp(:,1:2),'g*');
% plot(P2(:,1:2),'g*');
tform=affine2d(tform);
J=imwarp(w,tform);

figure;
imshow(J);


