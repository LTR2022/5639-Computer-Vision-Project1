%% parameter 

clc;
clear;          % clear all
 
tsigma=0.5;       %sigma of time temporal derivative gaussian filter
Threshold=0;    %Threshold of 1D time temporal derivative filter [-1 0 1]

SmoothfilterSize=5;    % size of smooth filter
ssigma=1;              % sigma of gaussian smoothing filter


%% load image
n=354;
I=cell(1,n+1);
for i=2:n
    imageName=strcat(['RedChair/RedChair/advbgst1_21_',sprintf('%04i',i)],'.jpg');
    I{i} = rgb2gray(imread(imageName));  
end
 
% boundry padding
I{1}= rgb2gray(imread(strcat('RedChair/RedChair/advbgst1_21_0002','.jpg')));
I{355}= rgb2gray(imread(strcat('RedChair/RedChair/advbgst1_21_0354','.jpg')));


%% main process

m=double(uint8(5*tsigma))+rem(double(uint8(5*tsigma))+1,2);
r1=uint16(m/2)-rem(m,2);
m=double(uint8(10*tsigma))+rem(double(uint8(10*tsigma))+1,2);
r2=uint16(m/2)-rem(m,2);
m=double(uint8(20*tsigma))+rem(double(uint8(20*tsigma))+1,2);
r3=uint16(m/2)-rem(m,2);

[DevImageBox,DevImageGaussian]= tfilterfunction(I,tsigma,Threshold);   % temporal derivative
[DevImageBox1,DevImageGaussian1]= tfilterfunction(I,2*tsigma,Threshold);
[DevImageBox2,DevImageGaussian2]= tfilterfunction(I,4*tsigma,Threshold);

SmoothImage1 = smooth_box(I,SmoothfilterSize);
[DevSmoothImage_Box1,DevSmoothImage_Gaussian1]= tfilterfunction(SmoothImage1,2*tsigma,Threshold); % temporal derivative again

SmoothImage2 = smooth(I,SmoothfilterSize,ssigma); % smooth all image one-by-one
[DevSmoothImage_Box2,DevSmoothImage_Gaussian2]= tfilterfunction(SmoothImage2,2*tsigma,Threshold); % temporal derivative again

% a
% for i= 2:350
%     subplot(3,2,1);
%     imshow(I{i});
%     title('raw image');
%     subplot(3,2,3);
%     imshow(DevImageBox{i});
%     title('simple 0.5[−1, 0, 1] filter');
%     subplot(3,2,4);
%     imshow(DevImageGaussian{i+r1-1});
%     title('gaussian filter with tsigma 0.5');
%     subplot(3,2,5);
%     imshow(DevImageGaussian1{i+r2-1});
%     title('gaussian filter with tsigma 1');
%     subplot(3,2,6);
%     imshow(DevImageGaussian2{i+r3-1});
%     title('gaussian filter with tsigma 2');
%     pause(0.1);
% end

% b box filter 
% need to change SmoothfilterSize
% for i= 2:350
%     subplot(3,2,1);
%     imshow(I{i});
%     title('raw image');
%     subplot(3,2,3);
%     imshow(DevImageBox{i});
%     title('simple 0.5[−1, 0, 1] filter');
%     subplot(3,2,4);
%     imshow(DevImageGaussian1{i+r2-1});
%     title('gaussian filter with tsigma 1');
%     subplot(3,2,5);
%     imshow(DevSmoothImage_Box1{i});
%     title('simple 0.5[−1, 0, 1], box filter');
%     subplot(3,2,6);
%     imshow(DevSmoothImage_Gaussian1{i+r2-1});
%     title('gaussian, tsigma 1, box filter');
%     pause(0.1);
% end

% b gaussian filter 
for i= 2:350
    subplot(3,2,1);
    imshow(I{i});
    title('raw image');
    subplot(3,2,3);
    imshow(DevImageBox{i});
    title('simple 0.5[−1, 0, 1] filter');
    subplot(3,2,4);
    imshow(DevImageGaussian1{i+r2-1});
    title('gaussian filter with tsigma 1');
    subplot(3,2,5);
    imshow(DevSmoothImage_Box2{i});
    title('simple 0.5[−1, 0, 1], gaussian smoothing');
    subplot(3,2,6);
    imshow(DevSmoothImage_Gaussian2{i+r2-1});
    title('gaussian, tsigma 1, gaussian smoothing');
    pause(0.1);
end
    

%% tfilter function
function [ DevImageBox, DevImageGaussian] = tfilterfunction(I,sigma,Threshold)
    
    % 1D box filter
    tfilterBox = 0.5*[-1, 0, 1];
    
    % gaussian filter
    m=double(uint8(5*sigma))+rem(double(uint8(5*sigma))+1,2);
    r=uint16(m/2)-rem(m,2);
    
    tfilterGaussian=fspecial('gaussian',[1,m],sigma);
    tfilterGaussian1 = tfilterGaussian;
    
    r = double(r);
    for j= -r:r   
            tfilterGaussian1(j+r+1)= -1*double((-j)/sigma^2)*double(tfilterGaussian(j+r+1));  
    end
    tfilterGaussian = tfilterGaussian1;
    Sum=sum(abs(tfilterGaussian));
    tfilterGaussian=1/Sum*tfilterGaussian;
    
    % output image initicizalize
    DevImageBox=cell(1,355);
    DevImageGaussian=cell(1,355+m);

    % padding
    for i=1:355
             DevImageBox{i} = zeros(240,320);
    end
    
    for i=1:355+2*r
            DevImageGaussian{i} = zeros(240,320);
    end
    
     % modify the size of input image for gaussian filter
    GaussianI=cell(1,355+m);
    for i=1:r
        GaussianI{i}= 0*I{i};
    end
    
    for i=r+1:355+r
        GaussianI{i}= I{i-r};
    end

    for i=355+r:355+2*r
        GaussianI{i}= 0*I{1};
    end    
    
    % filting with boxfilter
    for i= 2:354
         DevImageBox{i} = tfilterBox(3)*double(I{i+1})+tfilterBox(1)*double(I{i-1})+tfilterBox(2)*double(I{i});  
         DevImageBox{i} =abs(DevImageBox{i})-Threshold;
         %imshow(DevImageBox{i});
    end
    
    % r = uint8(r);    
    % filting with Gaussian filter
    for i =1+r:354+r
        for j= -r:r   
            DevImageGaussian{i}= DevImageGaussian{i}+ tfilterGaussian(j+r+1)*double(GaussianI{i+j});  
        end
        DevImageGaussian{i} =abs(DevImageGaussian{i})-Threshold;
        %imshow(DevImageGaussian{i});
        %imshow(uint8(DevImageGaussian{i}));
    end

end

%% smoothing function - gaussian
function [SmoothImage]= smooth(I,Size,ssigma)
    Smooth = fspecial('gaussian',Size,ssigma);
    for i=1:355
        SmoothImage{i}=imfilter(I{i},Smooth,'replicate');
    end
end

%% smoothing function - box
function [SmoothImage]= smooth_box(I,Size)
    Smooth = 1/(Size^2)*ones(Size,Size);
    for i=1:355
        SmoothImage{i}=imfilter(I{i},Smooth,'replicate');
    end
end





