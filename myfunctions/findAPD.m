function APD = findAPD(t_find,AP_find,APD_lim)

[~,locmaxAPD] = max(AP_find);
APD_cond = (max(AP_find) - min(AP_find)) * (1-APD_lim) + min(AP_find); 

if min(AP_find(1:locmaxAPD)) > APD_cond
    AP_start_loc = 1; 
else
    [~,AP_start_loc] = min(abs(AP_find(1:locmaxAPD)-APD_cond));
end

if  min(AP_find(locmaxAPD:end)) > APD_cond
    AP_end_loc = length(AP_find(locmaxAPD:end)); 
else
    [~,AP_end_loc] = min(abs(AP_find(locmaxAPD:end)-APD_cond));
end

APD = t_find(AP_end_loc+locmaxAPD) - t_find(AP_start_loc);
end

