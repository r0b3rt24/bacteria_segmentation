%% Pipeline
% Put the main procedure code here

files = dir(fullfile("data", '*.jpg'));

% files(1).name  % This is how you get the first image filename
ii = randperm(267,30)

for i = 1:30
    im = imread("./data/" + files(ii(i)).name);

    im = imgaussfilt(im,32);  % apply gauss filter to get rid of 
    
    imshow(im)
    im_after = im;

    for k = 0.176:0.16:0.816
        tem = im2bw(im, k);  % So the key here is to fin the level for all the images
        [w,l] = size(tem);
        if sum(tem(:)) < w*l*0.67 && sum(tem(:)) > w*l*0.33
            tem = bwareaopen(tem,round(w*l*0.2));
            tem = ~bwareaopen(~tem,round(w*l*0.001));
            if sum(sum(~bwareaopen(~tem,round(w*l*0.2)))) < w*l*0.33
                continue
            end
            im_after = tem;
            break
        end
    end

    [B,L]= bwboundaries(im_after);

    figure
    imshow(label2rgb(L, @jet, [.1 .1 .1]))
    hold on
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
    end
    [path, name] = fileparts(files(ii(i)).name);
    saveas(gcf,strcat(name,"_Seg.png"));
end

%% Single Run
L = bacteria_segment("./data/PIL-5_3dayLBCR-2.jpg");
B = edge(L);

%% Evaluation
L = bacteria_segment("./data/PIL-174_3dayLBCR-3.jpg");
B = edge(L);
ground = im2bw(imread("./GT/PIL-174_3dayLBCR-3_GT.jpg"),0.9);

B_g = edge(ground);

% calculating Hausdorff distance, we couldn't run this within 30mins. 

% dist = [];
% for i = 1:size(B, 1)
%     a = ones(size(B_g, 1), 1) * B(i, :);
%     b = (a - B_g) .* (a - B_g);
%     b = sqrt(b * ones(size(B,2),1));
%     dist(i) = min(b);
% endL
% 
% dist = max(dist);
% 
% Hausdorff = dist

overlap_area = length(find((-1* (ground-1) + L) == 2));
g_area = length(find(-1* (ground-1) == 1));
my_area = length(find(L == 1));

dice_coe = (2*overlap_area)/(g_area + my_area)

% below are code support visualization 

% figure;
% imshow(B);
% hold on;
% % plot(my_y(1),my_x(1), 'r*');
% axis on;
% figure;
% imshow(B_g);
% axis on;
% figure;
% imshow(B+B_g)

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

function L =bacteria_segment(data_path)
    im = imread(data_path);
    
    figure
    imshow(im)

    im = imgaussfilt(im,32);  % apply gauss filter to get rid of 
    
    im_after = im;

    for k = 0.176:0.16:0.816
        tem = im2bw(im, k);  % So the key here is to fin the level for all the images
        [w,l] = size(tem);
        if sum(tem(:)) < w*l*0.67 && sum(tem(:)) > w*l*0.33
            tem = bwareaopen(tem,round(w*l*0.2));
            tem = ~bwareaopen(~tem,round(w*l*0.001));
            if sum(sum(~bwareaopen(~tem,round(w*l*0.2)))) < w*l*0.33
                continue
            end
            im_after = tem;
            break
        end
    end

    [B,L]= bwboundaries(im_after, 'noholes');
    
    figure
    imshow(label2rgb(L, @jet, [.1 .1 .1]))
    hold on
    for k = 1:length(B)
       boundary = B{k};
       plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
    end
    
    % Uncomment this if you need to save the image.
%     [path, name] = fileparts(files(ii(i)).name);
%     saveas(gcf,strcat(name,"_Seg.png")); 
    
end