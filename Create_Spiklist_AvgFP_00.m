function [Spik_list, avgFP, n_avg] = Create_Spiklist_AvgFP_00(Header,stim_list,AP,AD)

un_stim   = unique(stim_list);
un_stim   = un_stim(find(~isnan(un_stim)));
n_stim    = length(un_stim);
stim_ind  = zeros(size(un_stim));

n_stim_rep = 0;
for i1 = 1:size(un_stim,1)
    i_stim = find(stim_list == un_stim(i1));
    n_stim_rep = max(max(size(i_stim)),n_stim_rep);
end

% CREATE A SPIKE LIST AND LFP AVG ACCORDING TO THE ONLINE FORMAT
% update the spik_list for plotting (1 stim, 2 n_rep, 3 ypos 
% ( (stim-1)+1/max(Header.nrep)*n_rep 4 chan 5 un 6 time
Spik_list = [];
dim_avgFP = size(AD);
avgFP     = zeros([dim_avgFP(1:end-1),n_stim]);

% CONVERT THIS FORMAT (AP)
% r: occurence; c: 1 - trial, 
% 2 - stimulus, 3 - ID-channel, 4 - ID-unit, 5 - ts 
% according to file start, 6 - ts according to stimulus-
% onset,
% INTO THIS FORMAT (Spik_list)
%(1 stim, 2 n_rep, 3
% ypos ( (stim-1)+1/max(Header.n_rep)*n_rep 4 chan 5 un 6
% time
if ~isempty(AP)
    Spik_list = zeros(size(AP,1),6);
    Spik_list(:,[1 4 5 6])  = AP(:,[2,3,4,6]); % fill up all the data possible without looping
end
for i1 = 1:size(stim_list,1)
    if ~isnan(stim_list(i1));
        stim_ind(find(stim_list(i1) == un_stim)) = stim_ind(find(stim_list(i1) == un_stim))+1;
        evt_id = find(un_stim == stim_list(i1));
        if ~isempty(AP)
            i_spik_trial = find(AP(:,1) == i1);
            n_spik_trial = size(i_spik_trial,1);
            if ~isempty(i_spik_trial) & i_spik_trial ~= 0
                try
                Spik_list(i_spik_trial,[2 3]) = repmat([stim_ind(find(stim_list(i1) == un_stim)), evt_id-1+(1/n_stim_rep)*(stim_ind(find(stim_list(i1) == un_stim))-1)], [n_spik_trial,1]);
                catch
                    keyboard;
                end
            end
        end
        if ~isempty(AD)
            try
                avgFP(:,:,evt_id) = (avgFP(:,:,evt_id)*((stim_ind(find(stim_list(i1) == un_stim))-1)/stim_ind(find(stim_list(i1) == un_stim))))+AD(:,:,i1)/stim_ind(find(stim_list(i1) == un_stim));  
            catch
                keyboard;
            end
        end
        clear n_spik_trial i_spik_trial evt_id
    end
end
n_avg = stim_ind;