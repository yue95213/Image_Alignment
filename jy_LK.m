%%%%%%%%%%% FA_LK %%%%%%%%%%%%%
clc;
clear;
close all;

input1=imread('tank.jpg');
input2=imread('tank3.jpg');
% template=rgb2gray(input1);
% image=rgb2gray(input2);
template=input1;
image=input2;

% imwrite(template,'00000001.jpg');
% imwrite(image,'0000002.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%initializing transformation matrix using corner features
% %%%%%%%%%% corners of template image %%%%%%%%%%%%%
% figure(1);
% imshow(template);
% title('template');
% C1=detectHarrisFeatures(image);
% hold on;
% plot(C1);
% 
% %%%%%%%%%% corners of warped image %%%%%%%%%%%%
% figure(2);
% imshow(image);
% title('image');
% C2=detectHarrisFeatures(image);
% hold on;
% plot(C2);

% %%%%%%%%%% matching features %%%%%%%%%%%
% %%Extract the neighborhood features.
% [features1,valid_points1] = extractFeatures(image,C1);
% [features2,valid_points2] = extractFeatures(image,C2);
% 
% %%Match the features.
% indexPairs = matchFeatures(features1,features2);
% 
% %%Retrieve the locations of the corresponding points for each image.
% matchedPoints1 = valid_points1(indexPairs(:,1),:);
% matchedPoints2 = valid_points2(indexPairs(:,2),:);
% 
% %%Get the affine transformation matrix of matched points between two images
% transform=estimateGeometricTransform(matchedPoints1,matchedPoints2,'affine');
% %%use 2*3 form of affine matrix
% tform=transform.T(1:2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%default parameters
par0.iterations =100;
par0.levels = 1;
par0.transform = 'homography';

if exist('par','var')
    if ~isstruct(par)
        error('iat_LucasKanade: the datatype of parameters is not a matlab struct');
    end
    
    if isfield(par,'initwarp') && ~isfield(par,'transform')
        error('iat_LucasKanade: when you initialize the warp, you should define the type of transform as well');
    end
    params = iat_merge_param(par0, par);
else
    params = par0;
end

if ~iat_is_transform(params.transform)
    error('iat_LucasKanade: unknown transform type. Check the field .transform in parameters structure');
end

if strcmpi(params.transform,'similarity')
    params.transform = 'affine';
    warning('iat_LucasKanade: Lukas-Kanade does not support similarity transform. Warp was automatically changed to affine')
end
    
transform = params.transform;    



if isfield(params,'initwarp')
    warp = params.initwarp;
    szw = size(warp);
    switch lower(transform)
        case 'translation'
            nop =2;
            if (szw(1)~=2 || szw(2)~=1)
                error('iat_LucasKanade: the warp matrix must be 2x1 for translation transform');
            end
        case 'euclidean'
            nop = 3;
            if (szw(1)~=2 || szw(2)~=3)
                error('iat_LucasKanade: the warp matrix must be 2x3 for euclidean transform');
            end
        case 'affine'
            nop=6;
            if (szw(1)~=2 || szw(2)~=3)
                error('iat_LucasKanade: the warp matrix must be 2x3 for affine transform');
            end
        case 'homography'
            nop = 8;
            if (szw(1)~=3 || szw(2)~=3)
                error('iat_LucasKanade: the warp matrix must be 3x3 for homography transform');
            end
    end
else
    switch lower(transform)
        case 'translation'
            warp = zeros(2,1);
            nop =2;
        case 'euclidean'
            nop = 3;
            warp = [eye(2) zeros(2,1)];
        case 'affine'
            nop=6;
            warp = [eye(2) zeros(2,1)];
%             warp=[-1 0 0;0 1 0];

        case 'homography'
            nop = 8;
            warp = eye(3);
    end
    
end

% warp=tform;
% nop=6;

levels = params.levels;
noi = params.iterations;

% Number of color channels for both image and template
sZi3 = size(image,3);
sZt3 = size(template,3);

% Color format validity check for image (RGB or gray-scale)
if sZi3>1
    if ((sZi3==2) || (sZi3>3))
        error('iat_LucasKanade: Unknown color image format: check the number of channels');
    else
        image=rgb2gray(uint8(image));
    end
end

% Color format validity check for image (RGB or gray-scale)
if sZt3>1
    if ((sZt3==2) || (sZt3>3))
        error('iat_LucasKanade: Unknown color image format: check the number of channels');
    else
        template = rgb2gray(uint8(template));
    end
end

% Converting template and image to doubles
template = double(template);
image = double(image);

%% pyramid images
% The following for-loop creates pyramid images in cells IM and TEMP with varying names
% The variables IM{1} and TEMP{1} are the images with the highest resoltuion

TEMP{1} = template;
IM{1} = image;

% Enable smoothing (optional)
% f = fspecial('gaussian',[7 7],.5);
% TEMP{1} = imfilter(template,f);
% IM{1} = imfilter(image,f);

for nol=2:levels
    IM{nol} = imresize(IM{nol-1},.5);
    TEMP{nol} = imresize(TEMP{nol-1},.5);
end

% in case of pyramid implementation, the initial transformation must be
% appropriately modified
for ii=1:levels-1
    warp=iat_warp_updown(warp, transform, 0);
end



if levels==1
    disp('Lucas-Kanade is running in single-level mode....');
else
    disp('Lucas-Kanade is running in multi-level mode....');
end

%% Run Lucas-Kanade algorithm for each level of pyramid
for nol=levels:-1:1   %%pyramid
    if levels>1
        fprintf('Level %d...', nol);
    end
    im = IM{nol};
    [vx,vy]=gradient(im);
    
    temp = TEMP{nol};
    
   [A,B]=size(temp);
    margin = 0; % no margin (enable a margin if you want)
    
    nx=margin+1:B-margin;
    ny=margin+1:A-margin;
    temp=double(temp(ny,nx,:));
   
    
    for i=1:noi    %%iteration
        
        %disp(['LucasKanade: Level: ' num2str(nol) ', Iteration: ' num2str(i)])
        %Image interpolation method
        str='linear'; % bilinear interpolation 
        %str='cubic'; % cubic ibterpolation
        
        wim = iat_inverse_warping(im, warp, transform, nx, ny, str); %inverse (backward) warping
        
        if (i == noi) % the algorithm is executed (noi-1) times
            break;
        end
        
        % Gradient Image interpolation (warped gradients)
        wvx = iat_inverse_warping(vx, warp, transform, nx, ny, str);
        wvy = iat_inverse_warping(vy, warp, transform, nx, ny, str);
        
        % Compute the jacobian of warp transform
        J = iat_warp_jacobian(nx, ny, warp, transform);
        
        % Compute the jacobian of warped image wrt parameters (steepest
        % descent image)
        
        G = iat_image_jacobian(wvx, wvy, J, nop);
        
        % Compute Hessian and its inverse
        C= G' * G;% C: Hessian matrix
        %i_C = inv(C);
        
        %show Hessian matrix
%         hh=zeros(300,300);
%         for i_h=1:6
%             for j_h=1:6
%                 if (i_h==5||i_h==6||j_h==5||j_h==6)
%                     hh(1+50*(i_h-1):50*i_h,1+50*(j_h-1):50*j_h)=C(i_h,j_h)*30;
%                 else
%                     hh(1+50*(i_h-1):50*i_h,1+50*(j_h-1):50*j_h)=C(i_h,j_h);
%                 end
%                 if(i_h==5||i_h==6)
%                     if(j_h==5||j_h==6)
%                         hh(1+50*(i_h-1):50*i_h,1+50*(j_h-1):50*j_h)=C(i_h,j_h)*500;
%                     end
%                 end
%             end
%         end
%         imshow(mat2gray(hh));
%         imwrite(mat2gray(hh),'show_hessian.jpg');
        
%         %%imshow hessian inverse
%         ihh=zeros(300,300);
%         tempC=inv(C);
%         tempC=mat2gray(tempC);
%         for i_h=1:6
%             for j_h=1:6
%                 ihh(1+50*(i_h-1):50*i_h,1+50*(j_h-1):50*j_h)=tempC(i_h,j_h)*3; 
%                 if(i_h==5||i_h==6)
%                     if(j_h==5||j_h==6)
%                         ihh(1+50*(i_h-1):50*i_h,1+50*(j_h-1):50*j_h)=tempC(i_h,j_h);
%                     end
%                 end
%             end
%         end
%         imshow(mat2gray(ihh));
        
        % Compute error vector
        imerror = temp - wim;
        
        % Compute the projection of error vector into Jacobian G
        Ge = G' * imerror(:);
        
        % Compute the optimum parameter correction vector
        delta_p = C\Ge;
        
        % Update parmaters
        warp = iat_warp_update(warp, delta_p, transform);
        
        
    end
    
    % modify the parameteres appropriately for next pyramid level
    if (nol>1)
        warp = iat_warp_updown(warp, transform,1);
    end
      if levels>1
    fprintf('Done\n');
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [wim, ones_map] = iat_inverse_warping(im, warp, transform, nx, ny, str); %inverse (backward) warping

figure;
imshow(uint8(template));
title('template');

figure;
imshow(uint8(image));
title('image');

figure;
imshow(mat2gray(wim));
title('warped image');

imwrite(mat2gray(wim),'warped image.png');