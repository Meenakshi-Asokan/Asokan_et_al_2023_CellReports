%Meenakshi M Asokan's code to compute stLFP using linear deconvolution
%Related to Figure 5 from Asokan et al., Cell Rep 2023
%September 2023

%Inputs - sound evoked rasters of units from one region and LFP from the
%other region (saved after notch filtering and common mode reference)

sta = cell(1,length(units));

parfor clu = 1:length(units)
    spike_rasters = units(clu).sound_raster;
    temp_sta = cell(num_chans,2);
    for ch = 1:num_chans
        lfp = cell(1,length(sounds));
        for sound = 1:2
            lfp{sound} = lfp_car{sound}(:,:,ch);% with common mode ref
        end
        %filtered version
        fs = 1000;
        LowCutoff = 0.5;
        HighCutoff = 300;
        [B,A] = butter(2,LowCutoff/(fs/2),'high');
        [C,D] = butter(2,HighCutoff/(fs/2),'low');
        trace_filt = cell(1,length(sounds));
        for i = 1:length(sounds)
            traces = lfp{i};
            filt_trace = [];
            for ii = 1:size(traces,1)
            filt_trace(ii,:) = filtfilt(B,A,traces(ii,:));
            filt_trace(ii,:) = filtfilt(C,D,filt_trace(ii,:));
            end
            trace_filt{i} = filt_trace;
        end
        
        sta_temp = cell(1,length(sounds));
        
        for sound = 1:length(sounds)
            raster = spike_rasters{sound};
            raster_sound_win = cat(2,zeros(size(raster,1),1000),raster(:,1000+(1:5000)),zeros(size(raster,1),1000));
            s = [];
            for row = 1:size(raster_sound_win,1)
                s = [s raster_sound_win(row,:)];
            end
            s = s';    
            trace_sound = trace_filt{sound};   
            trace_sound_sound_win = trace_sound(:,1:7000);
            y = [];
            for row = 1:size(trace_sound,1)
                y = [y trace_sound_sound_win(row,:)];
            end
            y = y';
            
            %Construct the design matrix
            n = 2000;
            m = length(s);
            X = zeros(m,n);%n - length of hdr h; m - length of the spikes s
            temp = s;
            for i= n/2:-1:1%1000 samples of hdr before the spike and 1000 samples after the spike
                temp = [temp(2:end);0];
                X(:,i) = temp;
            end
            temp = s;
            for i= n/2+1:n
                X(:,i) = temp;
                temp = [0;temp(1:end-1)];
            end
            
            % Make sparse
            X = sparse(X);
        
            % Normalize
            normfactor = sqrt(sum(X.^2)); % cant use norm because of sparsematrix. X = X./normfactor; %matlab 2016b and later
            X = bsxfun(@rdivide,X,normfactor);
        
            % Calculate
            Y = y;
            %units of beta - actual V (same as Y)    
            [beta,ISTOP,ITN] = lsmr(X,double(Y),[],10^-8,10^-8,[],400); % ISTOP = reason why algorithm has terminated, ITN = iterations
            stLFPest = beta;
            sta_temp{sound} = stLFPest;
        
        end
        temp_sta{ch,1}=sta_temp{cs_plus};
        temp_sta{ch,2}=sta_temp{cs_minus};    
    end
    sta{clu} = temp_sta;
    clu
end