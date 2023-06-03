%% parameter 

clc;
clear;          % clear all
 
tsigma=1;       %sigma of time temporal derivative gaussian filter
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
I{1}= rgb2gray(imread(strcat(['RedChair/RedChair/advbgst1_21_0002'],'.jpg')));
I{355}= rgb2gray(imread(strcat(['RedChair/RedChair/advbgst1_21_0354'],'.jpg')));


%% main process

m=double(uint8(5*tsigma))+rem(double(uint8(5*tsigma))+1,2);
r=uint16(m/2)-rem(m,2);

[DevImageBox,DevImageGaussian]= tfilterfunction(I,tsigma,Threshold);   % temporal derivative

newDevImageBox = DevImageBox;
% initalize mean and std value of each frame
mean__ = zeros(349);
std__ = zeros(349);
for i = 2:350
    mean_ = mean2(DevImageBox{i});
    mean__(i-1) = mean_;
    if mean_ < 5
        maxval = 10;
    else
        maxval = 1;
    end

    std_ = std2(DevImageBox{i});
    std__(i) = std_;
    newDevImageBox{i} = DevImageBox{i}-maxval*std_;

end

% plot mean and std value of each frame
% x=2:350;
% % z=mean__(x);
% z = std__(x);
% plot(x,z)
% xlabel( 'frame' );
% ylabel( 'std value' );

% plot function
for i= 2:350
    subplot(1,2,1);
    imshow(DevImageBox{i});
    title('box_filter');
    subplot(1,2,2);
    imshow(newDevImageBox{i+r-1});
    title('auto threshold');
    pause(0.1);
end
    

%% tfilter function
function [DevImageBox, DevImageGaussian] = tfilterfunction(I, sigma, Threshold)
    
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
    end
    
    % filting with Gaussian filter
    DevImageGaussian1=DevImageGaussian;
    DevImageGaussian2=DevImageGaussian;
    for i =1+r:354+r
        for j= -r:r   
            DevImageGaussian{i}= DevImageGaussian{i}+ tfilterGaussian(j+r+1)*double(GaussianI{i+j});  
        end
        DevImageGaussian{i} =abs(DevImageGaussian{i})-Threshold;
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





