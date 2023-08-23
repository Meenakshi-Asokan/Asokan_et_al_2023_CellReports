%Meenakshi M Asokan's code for creating a sample video demonstrating changes in facial motion energy
%Last Updated August 2023
%%
clear; clc; close all
% behav_vs_phys = 'P';
behav_vs_phys = 'B';

if strcmp(behav_vs_phys,'B')
    data_root = 'F:\Behavior\AAB';
else
    data_root = 'F:\EPHYS\AAP';
end
mouse_num = 6;
session_num = 1;
animalID = sprintf('AA%s%d',behav_vs_phys,mouse_num);
data_dir = fullfile(sprintf('//proteus/Polley_Lab/File Transfer/MMA/Rach Tosca/AA%s%d/Session %d',behav_vs_phys,mouse_num,session_num));
run_num = 1;
aviFolder = sprintf('//proteus/Polley_Lab/File Transfer/MMA/Rach pupil/AA%s%d/Session %d',behav_vs_phys,mouse_num,session_num);
cd(aviFolder);
% path = fullfile(data_root,sprintf('AA%s%d/DLC_processed_data/Session %d/Run%d',behav_vs_phys,mouse_num,session_num,run_num));
% mkdir(path);
% cd(path);
savingDir = fullfile(data_root,sprintf('AA%s%d/orofacial_data/Session %d/Run%d',behav_vs_phys,mouse_num,session_num,run_num));
mkdir(savingDir);
%%
file = dir('*.avi');
filename = file(1).name;
K=VideoReader(filename);
im1 = (read(K,10));
figure
imshow(im1(:,:,1));
% get the first ROI, which is posterior to the whisker pad
% here you need to preload the facial_ROI.mat to reuse the same ROIs
% use h = imrect(gca); to generate new ROIs and drag
% and place the rectangle
% h=imrect(gca, crop1);
h=imrect(gca);

% get the coordinates
Crop=getPosition(h);
crop_full(1,1)=Crop(1,1);
crop_full(1,2)=Crop(1,2);
crop_full(1,3)=Crop(1,3);
crop_full(1,4)=Crop(1,4);

%%
files = dir('*.avi');
roi_full = cell(1,length(files));
tic
%parallelize reading the video files 
% parfor i = 1:length(files)
parfor i = 1:2
% for i = 1
    filename = files(i).name;
    v = VideoReader(filename);
    numframes = v.NumberOfFrames;
    for ii = 1:numframes
        im = read(v,ii);
        im_gray = rgb2gray(im);
        roi_full{i}(:,:,ii) = imcrop(im_gray, crop_full);
    end 
    i
end
run_time = toc
% save('roi_fullface_raw.mat','roi_full'); 
%%
roi_face = [];
for i = 1:length(roi_full)
    roi_face = cat(3,roi_face,roi_full{i});
    i
end

move_face = abs(diff(double(roi_face), 1, 3));
%%
move_face = cat(3,zeros(size(move_face,1),size(move_face,2)),move_face);
framesinAvi = [];
for i = 1:2
    framesinAvi = [framesinAvi size(roi_full{i},3)];
end
cumframesAvi = [0 cumsum(framesinAvi)];
for i = 1:2
    move_face_split{i} = move_face(:,:,cumframesAvi(i)+(1:framesinAvi(i)));
end
%%
Facevideo.move_face_split = move_face_split;
Facevideo.roi.face = crop_full;
Facevideo.files = files;
cd(savingDir);
save([animalID,'_Facevideo.mat'], 'Facevideo')
%%
load('tl.mat');
%%
num_trials = length(tl.trials);
orofacial_trial_data = struct;
%so we basically have 1s baseline, 1s sound,10 s ITI, totally thats 12s*30 = 360 fps
for trial_num = 1:5
% for trial_num = 1:num_trials 
    start_st = 1;
    start_num_frames = length(tl.trials{trial_num}.states(start_st).frameInAVI);
    start_frames = (start_num_frames - 15 + 1):1:start_num_frames;%I want only 0.5 s (15 frames) preceding sound onset
    sound_st = 2;
    sound_num_frames = length(tl.trials{trial_num}.states(sound_st).frameInAVI);
    extra_sound_frames = sound_num_frames - 60; %I want only 2 s of sound (60 frames), sound starts after 0.5 s delay within this state
    end_st = 3;
    end_frames = 9.5*30 - extra_sound_frames; %9.5s;30fps - extra sound frames
    facial_movement = [];
    orofacial_trial_data(trial_num).sound = tl.trials{trial_num}.Group;
    if strcmp(orofacial_trial_data(trial_num).sound,'Noiseburst')
        orofacial_trial_data(trial_num).sound_index = 1;
        orofacial_trial_data(trial_num).sound_level = tl.trials{trial_num}.Sound.Noiseburst.Level.dB_SPL;
    elseif strcmp(orofacial_trial_data(trial_num).sound,'Ripple')
        orofacial_trial_data(trial_num).sound_index = 2;
        orofacial_trial_data(trial_num).sound_level = 80;       
    elseif strcmp(orofacial_trial_data(trial_num).sound,'up-ramp')
        orofacial_trial_data(trial_num).sound_index = 3;
        orofacial_trial_data(trial_num).sound_level = 80;    
    elseif strcmp(orofacial_trial_data(trial_num).sound,'down-ramp')
        orofacial_trial_data(trial_num).sound_index = 4;
        orofacial_trial_data(trial_num).sound_level = 80;    
    end
    
    frames = [tl.trials{trial_num}.states(start_st).frameInAVI(start_frames); tl.trials{trial_num}.states(sound_st).frameInAVI; tl.trials{trial_num}.states(end_st).frameInAVI(1:end_frames)];
    avi_num = [tl.trials{trial_num}.states(start_st).aviNum(start_frames); tl.trials{trial_num}.states(sound_st).aviNum; tl.trials{trial_num}.states(end_st).aviNum(1:end_frames)];
    [uni,ind] = unique(avi_num);
    for ii = 1:length(uni)
        if (ii == length(uni))
            index = length(frames);
        else
            index = ind(ii+1)-1;
        end
        facial_movement = [facial_movement Facevideo.move_face_split{uni(ii)+1}(:,:,frames(ind(ii):index))];
    end
    orofacial_trial_data(trial_num).facial_movement = facial_movement;

end
%%
%create videowriter with 30 fps
writerObj = VideoWriter('face_motion_exampletrial_80dB_noiseburst.avi');
writerObj.FrameRate = 30;
%open the videowriter
open(writerObj);
numframes_video = 300;%10 s video total
%write frames to the video
% for i = 1:size(orofacial_trial_data(1).facial_movement,3)
for i = 1:numframes_video
% for i = 1:2
    G = orofacial_trial_data(1).facial_movement(:,:,i);
%     imagesc(G);
    % Now make an RGB image that matches display from IMAGESC:
    C = colormap(hot);  % Get the figure's colormap.
    L = size(C,1);
    % Scale the matrix to the range of the map.
    Gs = round(interp1(linspace(min(G(:)),max(G(:)),L),1:L,G));
    H = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
%     imagesc(H);
    if (i == 30)
        speaker = dir('*.jpg');
        img = imresize(im2double(imread(speaker(1).name)),0.05);
        CC = colormap;  % Get the figure's colormap.
        LL = size(CC,1);
        % Scale the matrix to the range of the map.
        imgs = round(interp1(linspace(min(img(:)),max(img(:)),LL),1:LL,img));
        img_rgb = reshape(CC(imgs,:),[size(imgs) 3]); % Make RGB image from scaled.
        H_new = imfuse(H,img_rgb,'blend','Scaling','joint');
%         imagesc(H_new);
        writeVideo(writerObj,H_new);
    end
    writeVideo(writerObj,H);
end
close(writerObj);
%%
%plot frames

% for i = 1:size(orofacial_trial_data(1).facial_movement,3)
figure();
colormap(hot);
for i = 14+(1:60)
    G = orofacial_trial_data(1).facial_movement(:,:,i);
    h1 = subplot(6,10,i-14);
    imagesc(G);
    if (i == 30)
        C = colormap(hot);  % Get the figure's colormap.
        L = size(C,1);
        % Scale the matrix to the range of the map.
        Gs = round(interp1(linspace(min(G(:)),max(G(:)),L),1:L,G));
        H = reshape(C(Gs,:),[size(Gs) 3]); % Make RGB image from scaled.
        speaker = dir('*.jpg');
        img = imresize(im2double(imread(speaker(1).name)),0.1);
        CC = colormap;  % Get the figure's colormap.
        LL = size(CC,1);
        % Scale the matrix to the range of the map.
        imgs = round(interp1(linspace(min(img(:)),max(img(:)),LL),1:LL,img));
        img_rgb = reshape(CC(imgs,:),[size(imgs) 3]); % Make RGB image from scaled.
        H_new = imfuse(H,img_rgb,'blend','Scaling','joint');
        imagesc(H_new);
    end
    xticks([]);
    yticks([]);
    caxis([0 100]);
end
set(gcf,'color','w');
% export_fig(sprintf('AAB%d_Session%d_fullface_hotmouse_frames',mouse_num,session_num),'-png','-pdf')
%%
figure();
colormap(hot);
for i = 32
    G = orofacial_trial_data(1).facial_movement(:,:,i);
%     subplot(4,15,i-15);
    imagesc(G);
    xticks([]);
    yticks([]);
end
% colorbar();
set(gcf,'color','w');
% export_fig(sprintf('AAB%d_Session%d_fullface_hotmouse_frame32',mouse_num,session_num),'-png','-pdf')

%%
