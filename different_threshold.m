%% parameter 
clc;
clear;          % clear all
 
tsigma=0.5;       %sigma of time temporal derivative gaussian filter
Threshold=5;    %Threshold of 1D time temporal derivative filter [-1 0 1]

SmoothfilterSize=3;    % size of smooth filter
ssigma=1;              % sigma of gaussian smoothing filter

Swich= 'd';
%  a: do temporal derivative with origin image set
%  b: using box filter smooth the image and then do temporal derivative
%  c: using gaussian filter smooth the image and then do temporal derivative
%  d: try different threshold

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
r=uint16(m/2)-rem(m,2);

[DevImageBox,DevImageGaussian]= tfilterfunction(I,tsigma,Threshold);   % temporal derivative
[DevImageBox1,DevImageGaussian1]= tfilterfunction(I,2*tsigma,Threshold);
[DevImageBox2,DevImageGaussian2]= tfilterfunction(I,4*tsigma,Threshold);
imagewithRestult_DevImageBox= add_result(I, DevImageBox);     % add result
imagewithRestult_DevImageBox1 = add_result(I, DevImageBox1);
imagewithRestult_DevImageBox2 = add_result(I, DevImageBox2);
imagewithRestult_DevImageGaussian = add_result(I, DevImageGaussian);
imagewithRestult_DevImageGaussian1 = add_result(I, DevImageGaussian1);
imagewithRestult_DevImageGaussian2 = add_result(I, DevImageGaussian2);

SmoothImage1 = smooth_box(I,SmoothfilterSize);
[DevSmoothImage_Box1,DevSmoothImage_Gaussian1]= tfilterfunction(SmoothImage1,2*tsigma,Threshold); % temporal derivative again
imagewithRestult_box1= add_result(I, DevSmoothImage_Box1);
imagewithRestult_Gaussian1= add_result(I, DevSmoothImage_Gaussian1);

SmoothImage2 = smooth(I,SmoothfilterSize,ssigma); % gaussian ssigma = 1
[DevSmoothImage_Box2,DevSmoothImage_Gaussian2]= tfilterfunction(SmoothImage2,2*tsigma,Threshold); % threshold = 5
imagewithRestult_Box2= add_result(I, DevSmoothImage_Box2);
imagewithRestult_Gaussian2= add_result(I, DevSmoothImage_Gaussian2);

[DevSmoothImage_Box3,DevSmoothImage_Gaussian3]= tfilterfunction(SmoothImage2,2*tsigma,0.2*Threshold); % threshold = 1
imagewithRestult_Gaussian3= add_result(I, DevSmoothImage_Gaussian3);

[DevSmoothImage_Box4,DevSmoothImage_Gaussian4]= tfilterfunction(SmoothImage2,2*tsigma,5*Threshold); % threshold = 25
imagewithRestult_Gaussian4= add_result(I, DevSmoothImage_Gaussian4);

[DevSmoothImage_Box5,DevSmoothImage_Gaussian5]= tfilterfunction(SmoothImage2,2*tsigma,0); % no threshold
imagewithRestult_Gaussian5= add_result(I, DevSmoothImage_Gaussian5);

%% a plot
if Swich=='a'
    for i= 2:350
    gg1 = figure(1);
        set(gg1, 'Position', [100,100,1600,1200], 'Color', 'white');

        j=1;
        subplot(3,4,j);
        imshow(I{i});
        title('raw image');
        j=j+4; 
        
        subplot(3,4,j);
        imshow(DevImageBox{i});
        title('simple 0.5[-1, 0, 1] filter');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_DevImageBox{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevImageGaussian{i});
        title('gaussian filter with tsigma 0.5');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_DevImageGaussian{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevImageGaussian1{i});
        title('gaussian filter with tsigma 1');
        j=j+1;
         subplot(3,4,j);
        imshow(imagewithRestult_DevImageGaussian1{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevImageGaussian2{i});
        title('gaussian filter with tsigma 2');
        j=j+1;
         subplot(3,4,j);
        imshow(imagewithRestult_DevImageGaussian2{i});
        title('ADD Result');
        j=j+1;

        pause(1);
    end
end
%% b box filter 3*3 and 5*5
if Swich=='b'
    for i= 2:350
        gg1 = figure(1);
        set(gg1, 'Position', [100,100,1600,1200], 'Color', 'white');

        j=1;
        subplot(3,4,j);
        imshow(I{i});
        title('raw image');
        j=j+4; 

        subplot(3,4,j);
        imshow(DevImageBox{i});
        title('simple 0.5[-1, 0, 1] filter');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_DevImageBox{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevImageGaussian1{i});
        title('gaussian filter with tsigma 1');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_DevImageGaussian1{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevSmoothImage_Box1{i});
        title('simple 0.5[-1, 0, 1], 3*3 box');
        j=j+1;
         subplot(3,4,j);
        imshow(imagewithRestult_box1{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevSmoothImage_Gaussian1{i});
        title('gaussian, tsigma 1, 3*3 box');
        j=j+1;
         subplot(3,4,j);
        imshow(imagewithRestult_Gaussian1{i});
        title('ADD Result');
        j=j+1;

        pause(1);
    end
end
%% c gaussian filter ssigma = 1
if Swich=='c'
    for i= 2:350
        gg1 = figure(1);
        set(gg1, 'Position', [100,100,1600,1200], 'Color', 'white');

        j=1;
        subplot(3,4,j);
        imshow(I{i});
        title('raw image3');
        j=j+4; 

        subplot(3,4,j);
        imshow(DevImageBox{i});
        title('simple 0.5[-1, 0, 1] filter');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_DevImageBox{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevImageGaussian1{i});
        title('gaussian filter with tsigma 1');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_DevImageGaussian1{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevSmoothImage_Box2{i});
        title('simple 0.5[-1, 0, 1], gaussian smoothing');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_Box2{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevSmoothImage_Gaussian2{i});
        title('gaussian, tsigma 1, gaussian smoothing');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_Gaussian2{i});
        title('ADD Result');
        j=j+1;

        pause(1);
    end
end    

%% d different threshold
if Swich == 'd'
    for i= 2:350
        gg1 = figure(1);
        set(gg1, 'Position', [100,100,1600,1200], 'Color', 'white');

        j=1;
        subplot(3,4,j);
        imshow(I{i});
        title('raw image3');
        j=j+4;

        subplot(3,4,j);
        imshow(DevSmoothImage_Gaussian5{i});% no threshold
        title('no threshold');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_Gaussian5{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevSmoothImage_Gaussian3{i});% threshold = 1
        title('threshold = 1');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_Gaussian3{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevSmoothImage_Gaussian2{i});% threshold = 5
        title('threshold = 5');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_Gaussian2{i});
        title('ADD Result');
        j=j+1;

        subplot(3,4,j);
        imshow(DevSmoothImage_Gaussian4{i});% threshold = 25
        title('threshold = 25');
        j=j+1;
        subplot(3,4,j);
        imshow(imagewithRestult_Gaussian4{i});
        title('ADD Result');
        j=j+1;

        pause(0.5);
    end
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

    % make time frame same
    for i = 2:350
        DevImageGaussian{i} = DevImageGaussian{i+r-1};
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

%% add result to image
function [imagewithRestult]= add_result(I, result)
    for i=1:355

        imagewithRestult{i}=double(I{i})+double(result{i});
        imagewithRestult{i}=uint8(imagewithRestult{i});
    end
end







