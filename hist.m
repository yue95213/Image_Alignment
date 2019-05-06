clc;
clear;
close all;

temp=imread('tank.jpg');
temp=rgb2gray(temp);
im=imread('tank3.jpg');
im=rgb2gray(im);    
% figure;
% imshow(mat2gray(im));
% title('reference image');

%%pdf
[counts1,binLocations1] = imhist(temp);
sum1=sum(counts1);
pdf_temp=counts1/sum1;
[counts2,binLocations2] = imhist(im);
sum2=sum(counts2);
pdf_im=counts2/sum2;

%%intensity standardization
%template
intensity_sum1=counts1'*binLocations1;
mean1=intensity_sum1/sum1;
z1=binLocations1-mean1;
norm1=norm(double(temp(:)-mean1))/190;
z1=z1/norm1;
%image
intensity_sum2=counts2'*binLocations2;
mean2=intensity_sum2/sum2;
z2=binLocations2-mean2;
norm2=norm(double(im(:)-mean2))/190;
z2=z2/norm2;

%%image standardization
%template
s_temp=(double(temp)-mean1)./norm1;
%image
s_im=(double(im)-mean2)./norm2;

%%show graphs
fit_result1=fit(binLocations1,pdf_temp,'gauss1');
fit_result2=fit(binLocations2,pdf_im,'gauss1');
figure;
plot(fit_result1,'b');
hold on
plot(fit_result2,'r');
xlabel('pixel intensity');
ylabel('PDF');
% xlim([-0.005,0.005]);
axis([0 255 0 0.04]);
legend('Target','Input');

% %%target function
% a1 =     0.021 ;
% b1 =       105.5 ;
% c1 =       27.21 ;
% x1=0:0.1:255;
% func1=@(x)  a1*exp(-((x-b1)./c1).^2);
% figure;
% plot(x1,func1(x1));
% integral(func1,0,255)

%%s_target function
a1 =     0.01423 ;
b1 =     0;
c1 =   1.425  ;
func1=@(x)  a1*exp(-((x-b1)./c1).^2);
integral(func1,-Inf,Inf)
x1=-6:0.001:10;
figure;
plot(x1,func1(x1),'b');

%%input function
a1 =     0.01472 ;
b1 =      0;
c1 =       1.388  ;
func2 =@(x)  a1*exp(-((x-b1)./c1).^2);
x2=-6:0.001:10;
hold on;
plot(x2,func2(x2),'r');
xlabel('standardized pixel intensity');
ylabel('PDF');
legend('Target','Input');
