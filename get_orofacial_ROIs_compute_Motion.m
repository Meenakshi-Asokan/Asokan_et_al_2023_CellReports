%Meenakshi M Asokan's code for computing orofacial motion
%Last Updated August 2023
%%
clear; clc; close all
% behav_vs_phys = 'F';
behav_vs_phys = 'P';
% behav_vs_phys = 'B';
mouse_num = 23;
session_num = 1;
run_num = 1;
if strcmp(behav_vs_phys,'B')
    data_root = 'D:\Behavior\AAB';
    animalID = sprintf('AA%s%d',behav_vs_phys,mouse_num);

elseif strcmp(behav_vs_phys,'P')
    data_root = 'D:\EPHYS\AAP';
    animalID = sprintf('AA%s%d',behav_vs_phys,mouse_num);
else
    data_root = 'D:\FiberPhotometry\FAA';
    animalID = sprintf('%sAA%02d',behav_vs_phys,mouse_num);
    data_dir = fullfile(sprintf('//proteus/Polley_Lab/Meenakshi/DataTosca/%sAA%02d/Session %d',behav_vs_phys,mouse_num,session_num));
    aviFolder = sprintf('//proteus/Polley_Lab/Meenakshi/DataVideos/%sAA%02d/Session %d/Run%d',behav_vs_phys,mouse_num,session_num,run_num);
    savingDir = fullfile(data_root,sprintf('%sAA%02d/orofacial_data/Session %d/Run%d',behav_vs_phys,mouse_num,session_num,run_num));

end

cd(aviFolder);
mkdir(savingDir);
%%
file = dir('*.avi');
filename = file(1).name;
K=VideoReader(filename);
im1 = (read(K,10));
figure
imshow(im1(:,:,1));
% get the first ROI, which is posterior to the whisker pad
% here you can preload the facial_ROI.mat to reuse the same ROIs
% or use h = imrect(gca); to generate new ROIs and drag and place the rectangle
% h=imrect(gca, crop1);
% load('facial_ROI.mat');
h=imrect(gca);

% get the coordinates
Crop=getPosition(h);
crop1(1,1)=Crop(1,1);
crop1(1,2)=Crop(1,2);
crop1(1,3)=Crop(1,3);
crop1(1,4)=Crop(1,4);

% % get the whisker pad
% % h1=imrect(gca, crop2);
% h1=imrect(gca);
% 
% % get the coordinates
% Crop_new=getPosition(h1);
% crop2(1,1)=Crop_new(1,1);
% crop2(1,2)=Crop_new(1,2);
% crop2(1,3)=Crop_new(1,3);
% crop2(1,4)=Crop_new(1,4);
%%
%can also get pinna, body(?)
%% save the data if running for the first time
cd(savingDir);
save('facial_ROI.mat','crop1'); 
% save('facial_ROI_fullface.mat','crop1'); 
figure
im1_gray = rgb2gray(im1);
imshow(im1_gray);
hold on                                                                                 
rectangle('Position',crop1,'EdgeColor','b','LineWidth',2,'LineStyle','--'); hold on; % add the roi
export_fig('ROI',  '-png')
 
%%
cd(aviFolder);
files = dir('*.avi');
roi = cell(1,length(files));
tic
%parallelize reading the video files 
% for i = 1:length(files)
parfor i = 1:length(files)
    filename = files(i).name;
    v = VideoReader(filename);
    numframes = v.NumberOfFrames;
    for ii = 1:numframes
        im = read(v,ii);
        im_gray = rgb2gray(im);
        roi{i}(:,:,ii) = imcrop(im_gray, crop1);
%         roi_anterior{i}(:,:,ii) = imcrop(im_gray, crop2);
    end 
    i
end
toc
% save('roi_face_raw.mat','roi'); 
%%
% roi_face = [];
% for i = 1:length(roi)
%     roi_face = cat(3,roi_face,roi{i});
%     i
% end
% 
% move_face = abs(diff(double(roi_face), 1, 3));
% Mtrace_face = mean(squeeze(mean(move_face,1)),1);

%%
%compute diff between frames without concatenating everything to save
%computational time
tic
mean_diff_face = cell(1,length(roi));
for i = 1:length(roi)
diff_face = abs(diff(double(roi{i}), 1, 3));
if (i < length(roi))
left = double(roi{i}(:,:,end));
right = double(roi{i+1}(:,:,1));
inbetween_diff = abs(right - left);
diff_face = cat(3,diff_face,inbetween_diff);
end
mean_diff_face{i} = mean(squeeze(mean(diff_face,1)),1);
i
end
toc
Mtrace_face = [];
for i = 1:length(roi)
    Mtrace_face = [Mtrace_face,mean_diff_face{i}];
end
%%
Mtrace_face = [0 Mtrace_face];
Mtrace_face_split = cell(1,length(files));
framesinAvi = [];
for i = 1:length(files)
    framesinAvi = [framesinAvi size(roi{i},3)];
end
cumframesAvi = [0 cumsum(framesinAvi)];
for i = 1:length(files)
    Mtrace_face_split{i} = Mtrace_face(cumframesAvi(i)+(1:framesinAvi(i)));
end



% % %if appending to an already saved variable:
% cd(savingDir);
% load([animalID,'_Orofacial.mat']);

Orofacial.Mtrace_face = Mtrace_face;
Orofacial.Mtrace_face_split = Mtrace_face_split;
Orofacial.roi.face = crop1;


cd(savingDir);
save([animalID,'_Orofacial.mat'], 'Orofacial')