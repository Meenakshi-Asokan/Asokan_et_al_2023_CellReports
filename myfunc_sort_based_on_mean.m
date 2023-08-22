%display after sorting - template
%input - data: eg. psth for all units
%output - psth with units sorted based on average value of psth
function sorted_data = myfunc_sort_based_on_mean(data)
    avg_resp = [];
    for clu = 1:size(data,1)
        psth = data(clu,:);
        avg_resp(clu) = mean(psth);
    end

    [out,idx] = sort(avg_resp,'descend');
    sorted_data = [];
    for i = 1:size(data,1)
        sorted_data(i,:) = data(idx(i),:);
    end
end