clear all
clc
close all

% Extraction of data from a graph image
cd './'
image = ['Craig Cox 1.jpg'];
% load image and select axes rectangle
Img = imread(image);
imshow(Img,'InitialMagnification','fit')
title('Select graph area from origin to axis maxima!','fontsize',14);
set(gcf, 'Position', get(0,'Screensize'));                % Maximise figure
pause(1)
frame = getrect;
X1 = round(frame(1));
Y1 = round(frame(2));
dX = ceil(frame(3)/2)*2;  % even X width
dY = ceil(frame(4)/2)*2;  % even Y width
X2 = X1 + dX - 1;
Y2 = Y1 + dY - 1;

% Get absolute axis values
Xmin = input('Enter X origin: ');
Xmax = input('Enter X maximum: ');
Ymin = input('Enter Y origin: ');
Ymax = input('Enter Y maximum: ');

figure(1)
title('Klick on datapoint or outside graph to exit','fontsize',14)

exit=0; i=0; Data = []; j = 1;                       % initialize variables
while exit==0
 i=i+1;
 [x y] = ginput(1);
 if x<X2 && y<Y2 && x>X1 && y>Y1
     Data{j}(i,:) = [x y];
 else
     close
     line = input('Additional set of datapoints (y/n)? \n','s');
     if line=='y'
         j = j+1;
         i=0;
         imshow(Img,'InitialMagnification','fit')
         title('Klick on datapoint or outside graph to exit','fontsize',14)
         set(gcf, 'Position', get(0,'Screensize'));       % Maximize figure
     else
         exit=1;
     end
 end
end

% Dimentionalize Datapoints

for m=1:j
    % pixel coordinates
    xpix = Data{m}(:,1);
    ypix = Data{m}(:,2);
    % true coordinates
    ov = ones(length(xpix),1);
    xtrue = ov*Xmin+(xpix-X1*ov)/dX*(Xmax-Xmin);
    ytrue = ov*Ymin+(dY*ov-ypix+ov*Y1)/dY*(Ymax-Ymin);
    Data{m} = [xtrue ytrue];
    figure(2)
    plot(xtrue,ytrue,'rx')
    hold on
end
grid
hold off