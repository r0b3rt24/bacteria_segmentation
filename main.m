%% Pipeline
% Put the main procedure code here

files = dir(fullfile("data", '*.jpg'));
% TODO: 1. find the threshold 2. find the largest circle in B

% files(1).name  % This is how you get the first image filename
im = imread("./data/" + files(67).name);
im = imgaussfilt(im,32);
% im = im2grey(im);

% Find the threshold

im = im2bw(im, 0.5);  % So the key here is to fin the level for all the images

[B,L]= bwboundaries(im);
% [x,y] = size(im)
% im = regiongrowing(im,floor(x/2),floor(y/2))
% [centers, radii] = imfindcircles(im, [3000,6000], 'ObjectPolarity', 'dark')

imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

%% Evaluation
ground = im2bw(imread("ground_truth.png"));
B = edge(seg);
orig_size = size(B);

ground = imresize(ground, orig_size, 'nearest');

B_g = edge(im2bw(ground));

dist = [];
for i = 1:size(B, 1)
    a = ones(size(B_g, 1), 1) * B(i, :);
    b = (a - B_g) .* (a - B_g);
    b = sqrt(b * ones(size(B,2),1));
    dist(i) = min(b);
end

dist = max(dist);

Hausdorff = dist

overlap_area = length(find((-1* (ground-1) + seg) == 2));
g_area = length(find(-1* (ground-1) == 1));
my_area = length(find(seg == 1));

dice_coe = (2*overlap_area)/(g_area + my_area)

figure;
imshow(B);
hold on;
% plot(my_y(1),my_x(1), 'r*');
axis on;
figure;
imshow(B_g);
axis on;
figure;
imshow(B+B_g)

% figure;
% imshow(im)

% figure;
% im = imread("./data/" + files(33).name);
% im = im2bw(im, 0.5);  % So the key here is to fin the level for all the images
% 
% imshow(im)



%% Functions
% Put the suppoert function code here

function J=regiongrowing(I,x,y,reg_maxdist)
    if(exist('reg_maxdist','var')==0), reg_maxdist=0.2; end
    if(exist('y','var')==0), figure, imshow(I,[]); [y,x]=getpts; y=round(y(1)); x=round(x(1)); end
    J = zeros(size(I));
    Isizes = size(I);
    reg_mean = I(x,y);
    reg_size = 1;
    neg_free = 10000; neg_pos=0;
    neg_list = zeros(neg_free,3); 
    pixdist=0;
    neigb=[-1 0; 1 0; 0 -1;0 1];
    while(pixdist<reg_maxdist&&reg_size<numel(I))
        for j=1:4,
            xn = x +neigb(j,1); yn = y +neigb(j,2);
            ins=(xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2));
            if(ins&&(J(xn,yn)==0)) 
                    neg_pos = neg_pos+1;
                    neg_list(neg_pos,:) = [xn yn I(xn,yn)]; J(xn,yn)=1;
            end
        end
        if(neg_pos+10>neg_free), neg_free=neg_free+10000; neg_list((neg_pos+1):neg_free,:)=0; end
        dist = abs(neg_list(1:neg_pos,3)-reg_mean);
        [pixdist, index] = min(dist);
        J(x,y)=2; reg_size=reg_size+1;

        reg_mean= (reg_mean*reg_size + neg_list(index,3))/(reg_size+1);

        x = neg_list(index,1); y = neg_list(index,2);

        neg_list(index,:)=neg_list(neg_pos,:); neg_pos=neg_pos-1;
    end
    J=J>1;
end