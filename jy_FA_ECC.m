%%%%%%%%%%% FA_ECC %%%%%%%%%%%%%
clc;
clear;
close all;

input1=imread('tank.jpg');
input2=imread('tank3.jpg');
template=rgb2gray(input1);
image=rgb2gray(input2);
% template=input1;
% image=input2;

figure(1);
imshow(template);
title('template');

figure(2);
imshow(image);
title('image');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [warp,rho]=jy_alignment(image.temp,par1);
%default parameters
par0.iterations =50;
par0.levels = 1;
par0.transform = 'affine';

%%%%%% rho in every iteration %%%%%%
rho_1=zeros(par0.iterations,1);

if exist('par','var')  %%check existence of 'par',check only for variables. (return 0"_not exist" / 1 "_variables")
    if ~isstruct(par)  %%returns 1 or 0 if par is a MATLAB structure or not.
        error('iat_ecc: the datatype of parameters is not a matlab struct');
    end
    
    if isfield(par,'initwarp') && ~isfield(par,'transform') %%check whether 'initwarp'/'transform' is the field of structure 'par'
        error('iat_ecc: when you initialize the warp, you should define the type of transform as well');
    end
    params = iat_merge_param(par0, par);  %%merge parameter sets 'par0' and 'par' into new parameters 'params'
else
    params = par0;
end

if ~iat_is_transform(params.transform) %%check if a valid transformation name
    error('iat_ecc: unknown transform type. Check the field .transform in parameters structure');
end

if strcmpi(params.transform,'similarity')  %%compare strings ignoring case (same/different return 1/0)
    params.transform = 'affine';
    warning('iat_ecc: ECC does not support similarity transform. Warp was automatically changed to affine')
end
    
transform = params.transform;    
    
if isfield(params,'initwarp')
    warp = params.initwarp;
    szw = size(warp);
    switch lower(transform) %%convert string to lower case (see definition about transform above)
        case 'translation'
            nop =2;
            if (szw(1)~=2 || szw(2)~=1)
                error('iat_ecc: the warp matrix must be 2x1 for translation transform');
            end
        case 'euclidean'
            nop = 3;
            if (szw(1)~=2 || szw(2)~=3)
                error('iat_ecc: the warp matrix must be 2x3 for euclidean transform');
            end
        case 'affine'
            nop=6;
            if (szw(1)~=2 || szw(2)~=3)
                error('iat_ecc: the warp matrix must be 2x3 for affine transform');
            end
        case 'homography'
            
            nop = 8;
            if (szw(1)~=3 || szw(2)~=3)
                error('iat_ecc: the warp matrix must be 3x3 for homography transform');
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
        case 'homography'
            nop = 8;
            warp = eye(3);
    end
    
end


break_flag=0;

levels = params.levels;
noi = params.iterations;

% Number of color channels for both image and template
sZi3 = size(image,3);  %return the length of dimension '3', if image is rgb format, returns 3
sZt3 = size(template,3);


% Color format validity check for image (RGB or gray-scale)
if sZi3>1
    if ((sZi3==2) || (sZi3>3))
        error('iat_ecc: Unknown color image format: check the number of channels');
    else
        image=rgb2gray(uint8(image));
    end
end

% Color format validity check for image (RGB or gray-scale)
if sZt3>1
    if ((sZt3==2) || (sZt3>3))
        error('iat_ecc: Unknown color image format: check the number of channels');
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

for nol=2:levels  %%down-sampling
    IM{nol} = imresize(IM{nol-1},.5);
    TEMP{nol} = imresize(TEMP{nol-1},.5);
end

% in case of pyramid implementation, the initial transformation must be
% appropriately modified
for ii=1:levels-1
    warp=iat_warp_updown(warp, transform, 0);   %%'0/1' represents change transformation to adapt low/high resolution
end

if levels==1
    disp('ECC is running in single-level mode....');
else
    disp('ECC is running in multi-level mode....');
end

%% Run ECC algorithm for each level of pyramid
for nol=levels:-1:1  %%from low resolution to high resolution
    if levels>1
        fprintf('Level %d...', nol); %%start of this pyramid level
    end
    
    im = IM{nol};
    [vx,vy]=gradient(im);
    
    temp = TEMP{nol};
    
    [A,B]=size(temp);
    % Warning for tiny images
    if prod([A,B])<400  %%product of the elements of the matrix [A,B]
        disp(' -> ECC Warning: The size of images in high pyramid levels is quite small and it may cause errors.');
        disp(' -> Try fewer levels or larger images to avoid such errors.');
        disp(' -> Press any key to continue.') %%diplay the array
        pause
    end
    
    margin = 0; % no margin
    
    % Uncomment the follwing lines if you want to consider margins
    % m0 = mean([A,B]);
    % margin = floor(m0*.05/(2^(nol-1))); % 5-percent of the mean of [height,width]
    
    nx=margin+1:B-margin;
    ny=margin+1:A-margin;
    
%     nx=margin+150:B-150;
%     ny=margin+150:A-150;
    temp=double(temp(ny,nx,:)); %get new template considering margin
    
    % Compute the jacobian of warp transform (2*N???)
    J = iat_warp_jacobian(nx, ny, warp, transform);    %%invariant in affine case
    
    
    %% ECC, Forwards Additive Algorithm -------------------------------
    for i=1:noi  %%iteration
        
        %disp(['ECC: Level: ' num2str(nol) ', Iteration: ' num2str(i)])
        %num2str:convert numbers to character array 
        
        %Image interpolation method
        str='linear'; % bilinear interpolation
        %str='cubic'; % bicubic interpolation
        
        [wim, ones_map] = iat_inverse_warping(im, warp, transform, nx, ny, str); %inverse (backward) warping
        %%%%%%%%%%%%%interpolation??????????
        %get warped interpolated image 'wim' and binary double image
        %'ones_map' (areas in 'wim' come from supported area(1) or outside the borders(0))
        
        % consider the overlap for the zero-mean process
        numOfElem = sum(ones_map(:)~=0);  %sum number of pixels(/elements) that come from supported area 
        meanOfWim = sum(wim(ones_map~=0))/numOfElem; 
        %get mean of these elements above in 'wim'
        meanOfTemp = sum(temp(ones_map~=0))/numOfElem;  %... in 'template'
        
        % Compute zero-mean images
        wim = wim-meanOfWim;% zero-mean image; is useful for brightness change compensation, otherwise you can comment this line
        tempzm = temp-meanOfTemp; % zero-mean template
        
        wim(ones_map==0) = 0; % reject pixels outside the overlap
        tempzm(ones_map==0) = 0;
        
        normOfwim = norm(wim(:));  %return 2-norm of 'wim'
        % Save current correlation coefficient
        rho = dot(temp(:),wim(:)) / norm(tempzm(:)) / normOfwim;
        
        rho_1(i)=rho;
        
        if (i == noi) % the algorithm is executed (noi-1) times (/now is the 'noi' time iteration) 
            break;
        end
        
        % Gradient Image interpolation (warped gradients)  (K*2)
%         wvx = iat_inverse_warping(vx, warp, transform, nx, ny, str);
%         wvy = iat_inverse_warping(vy, warp, transform, nx, ny, str);
        [wvx,wvy]=gradient(wim);
        
%         % Compute the jacobian of warp transform (2*N???)
%         J = iat_warp_jacobian(nx, ny, warp, transform);    %%invariant in affine case
        
        % Compute the jacobian of warped image wrt parameters (matrix G in
        % the paper)  (K*N)
        G = iat_image_jacobian(wvx, wvy, J, nop);
        
        % Compute Hessian and its inverse
        C= G' * G;% C: Hessian matrix
        
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
        
         %%imshow hessian inverse
%         ihh=zeros(300,300);
%         tempC=inv(C);
%         tempC=mat2gray(tempC);
%         for i_h=1:6
%             for j_h=1:6
%                 ihh(1+50*(i_h-1):50*i_h,1+50*(j_h-1):50*j_h)=tempC(i_h,j_h)*7; 
%                 if(i_h==5||i_h==6)
%                     if(j_h==5||j_h==6)
%                         ihh(1+50*(i_h-1):50*i_h,1+50*(j_h-1):50*j_h)=tempC(i_h,j_h);
%                     end
%                 end
%             end
%         end
%         imshow(mat2gray(ihh));
        
        
        con=cond(C);  %condition number
        if con>1.0e+15
            disp('->ECC Warning: Badly conditioned Hessian matrix. Check the initialization or the overlap of images.')
        end
        %i_C = inv(C);
        
        % Compute projections of images into G
        Gt = G' * tempzm(:);
        Gw = G' * wim(:);
        
        
        %% ECC closed form solution
        
        vector = C\Gw; % inv(C)*Gw
        
        % Compute lamda parameter
        num = (normOfwim^2 - dot(Gw,vector));
        den = (dot(tempzm(:),wim(:)) - dot(Gt,vector));
        lambda = num / den;
        
        % Compute error vector
        imerror = lambda * tempzm - wim;
        
        % Compute the projection of error vector into Jacobian G
        Ge = G' * imerror(:);
        
        % Compute the optimum parameter correction vector
        delta_p = C\Ge; % inv(C)*Ge
        
        if (sum(isnan(delta_p)))>0 %Hessian is close to singular
            disp([' -> ECC algorithms stopped at ' num2str(i) '-th iteration of ' num2str(nol) '-th level due to bad condition of Hessian matrix.']);
            break_flag=1;
            break;  %%break from this iteration
        end
        
        % Update parmaters
        warp = iat_warp_update(warp, delta_p, transform);
        
        
    end   %%end of iteration
    
    if break_flag==1
        break;  %%break from this pyramid level
    end
    
    % modify the parameteres appropriately for next pyramid level
    if (nol>1)&&(break_flag==0)
        warp = iat_warp_updown(warp, transform,1);
    end
    if levels>1
    fprintf('Done\n');  %%end of this pyramid level
    end
end


if break_flag==1 % this conditional part is only executed when algorithm stops due to Hessian singularity
    for jj=1:nol-1
        warp = iat_warp_updown(warp, transform,1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[wim, ones_map] = iat_inverse_warping(im, warp, transform, nx, ny, str); %inverse (backward) warping
figure(3);
imshow(mat2gray(wim));
title('final warped image');
imwrite(mat2gray(wim),'ecc_warp.jpg');

%%%%%%%%%%% plot rho %%%%%%%%%%
figure(4);
plot(rho_1);