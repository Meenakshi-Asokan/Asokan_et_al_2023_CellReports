%Meenakshi M Asokan and Ke Chen's code
%Last updated August 2023
function data = DLC_Video_PreProcess_AnteriorPosteriorPupillabelsOnly(path)
cd(path)
data = struct
file = dir('*.csv');
for i = 1:length(file)
    T = readtable(file(i).name);
    data(i).name=file(i).name;
    data(i).Variables = T([1,2],:);
    temp = load_DLC_csv(file(i).name);
    
    % the x position: 2:3:end; y position: 3:3:end; likelihood position:
    % 4:3:end
    
    for j = 1:size(temp,1) % loops through each frame
        loc(j).x = temp(j, 2:3:end);
        loc(j).y = temp(j, 3:3:end);
        loc(j).likelihood = temp(j, 4:3:end);
        
    end
    
    % separate the pupil and snout
    for j = 1:size(temp, 1)
        data(i).pupil(j).x = loc(j).x(1: end-1);
        data(i).pupil(j).y = loc(j).y(1: end-1);
        data(i).pupil(j).likelihood = loc(j).likelihood(1: end-1);
    end
    clear loc
    
end
%%
% get the distance based on anterior and posterior marker of pupil
% anterior: is the max x location and posterior is the min x location
for i = 1:length(data)
    fprintf('Processing Video # %d\n', i)
    pupil = data(i).pupil; 
    low_likelihood_frames = [];
    for j = 1:length(pupil)
        pupil_data = pupil(j);    % loop through each frame
        [~, ind_anterior] = max(pupil_data.x);
        [~, ind_posterior] = min(pupil_data.x);
        % get the pupil 'diameter' using the anterior, posterior markers
        distance_ap(j) = sqrt((pupil_data.x(ind_anterior) - pupil_data.x(ind_posterior)).^2 + ...
            (pupil_data.y(ind_anterior) - pupil_data.y(ind_posterior)).^2);

        cutoff = 0.7; % frame with marker likelihood bigger than 0.7 was used, others will be interpolated
        indx = find(pupil_data.likelihood<cutoff);
        if ~(isempty(indx))
            low_likelihood_frames = [low_likelihood_frames j];
            fprintf('Frame # %d has likelihood less than %1.1f \n', j, cutoff)
        end
        
    end
    distance_ap_likelihood_corrected = distance_ap;
    for ind = 1:length(low_likelihood_frames)
        index = low_likelihood_frames(ind);
        lower_range = 1:index-1;            
        upper_range = index+1:length(pupil);
        lower_lim = max(setdiff(lower_range, low_likelihood_frames));
        upper_lim = min(setdiff(upper_range, low_likelihood_frames));
        if (isempty(lower_lim)&&isempty(lower_lim))
            distance_ap_likelihood_corrected(index) = distance_ap(index);%no way to correct(?) - trial should not be used(?)
        elseif isempty(lower_lim)
            distance_ap_likelihood_corrected(index) = distance_ap(upper_lim);
        elseif isempty(upper_lim)
            distance_ap_likelihood_corrected(index) = distance_ap(lower_lim);
        else
            distance_ap_likelihood_corrected(index) = (distance_ap(lower_lim)+distance_ap(upper_lim))/2;
        end
    end
    data(i).distance_ap = distance_ap;
    data(i).distance_ap_likelihood_corrected = distance_ap_likelihood_corrected;
end
% for i = 1:length(data)
%     fprintf('Processing Video # %d\n', i)
%     [distance_ap, pupil] = DLC_pupil(data(i).pupil);
%     data(i).pupil_fit = pupil;
%     data(i).distance_ap = distance_ap;
%     clear distance_ap pupil
% end
save('Summary_data', 'data')