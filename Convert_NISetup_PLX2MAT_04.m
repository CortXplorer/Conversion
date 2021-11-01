function  [AD,AP,stim_list] = Convert_NISetup_PLX2MAT_04(filename,save_choi,type_choi,t_pre,dur);
% Convertion here safes values in V!

% check input
if ~exist('filename') | strmatch(filename,'')
    [n,p] = uigetfile('*.plx');
else
    p = filename(1:max(strfind(filename,'\')));
    n = filename(max(strfind(filename,'\'))+1:end);
end

if ~exist('save_choi') | isempty(save_choi)
    save_choi = 1;
end

if ~exist('type_choi') | isempty(type_choi)
    type_choi = 1;
end

if ~exist('t_pre') | isempty(t_pre)
    t_pre   = .2;
    dur = .3+t_pre; % duration of trial including pre and post time
end

% stop execution if no file has been chosen
if n == 0; return; end

% inform about file name
fprintf('\n\n');
disp('*************************************************');
fprintf('PROCESSING          : %s\n',n);
fprintf('IN FOLDER           : %s\n\n',p);

% Calc file size of PLX file
f = dir(p);
size_file = f(strmatch(n,char(f.name))).bytes;
fprintf('SIZE OF FILE          : %.1f KB, %.2f MB\n', size_file/1024, size_file/1024/1024);

% read general data from the plx-file
% general
    [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = mexPlex(13,[p,n]);
    % get information about number of units, number of spikes, number of samples on AD chan ...
    [tscounts, wfcounts, evcounts, contcounts ] = mexPlex(4,[p,n], 1); % get info about enabled channels and units, newest mexw also gives info of continuous chans
% AD-Data
    [n_chan_tot,freqs]   = mexPlex(12,[p,n]); % output is sampling frequency; n_chan is number of channels the system supports
    [n_chan_tot,gains_ad]= mexPlex(11,[p,n]); % output is total gain (including preamp gain and digital gain)
    [n_chan_tot,names]   = mexPlex(15,[p,n]);
    try
        % plexon now has a way of checking how many sampling points have
        % been acquired per channel thus if n_sp > 0 that channel was
        % active; this works starting with mexPlex version 1.5.0.0
        [n_chan_tot, n_sp] = mexPlex(23,[p,n]); % plx_adchan_samplecounts
        ad_chan_on = find(n_sp ~= 0);
        Fs = freqs(ad_chan_on);
        n_ad_chan_on = length(ad_chan_on);
    catch
        % way to check which channels had been enabled for older mexPlex
        preampgain = 1000;
        ad_chan_on = find(gains_ad > preampgain);
        % alternative way: test which channels have more then 0 samples
        ad_chan_on  = find(evcounts > 0 & [1:length(evcounts)] >= 300 &    [1:length(evcounts)] <= 299+n_chan_tot)-299;     
        n_ad_chan_on = length(ad_chan_on);
        Fs = freqs(1);
    end
% AP-Data
    [n_chan_tot,thresh]  = mexPlex(9,[p,n]);   % output is array of threshold vals for each channel; plx_chan_thresholds
    [n_chan_tot,gains_ap]= mexPlex(8,[p,n]);    % output is total gain (including preamp gain and digital gain);plx_chan_gains
    [n_chan_tot,names]   = mexPlex(14,[p,n]); % plx_chan_names
    [n_chan_tot,filters] = mexPlex(10,[p,n]); % check which spikechannels have been filtered;plx_chan_filters
    % preliminary way to check which channels had been enabled (write email
    % to plexonincs
    mingain = 1000;
    ap_chan_on = find(gains_ap > mingain); n_ap_chan_on = length(ap_chan_on);
    

    [dummy_row, dummy_col] = find(tscounts);
    ID = [dummy_col-1, dummy_row];
    n_un_tot = size(ID,1);
       
% TRIGGER events in file
try
    [nev, evts, sv] = mexPlex(3,[p,n], 2); 

catch
    nev = 0; evts= 0; sv = 0;
end
n_trig = nev;

ts_trig = evts;

% STROBE events in file
spontaneous     = 1;
try
    [nev, evts, sv] = mexPlex(3,[p,n], 257);
    spontaneous     = 0;
    if sv == -1 | all(sv == 0) % there are no strobes, must be spont
        spontaneous     = 1;
    end
catch
    spontaneous     = 1;
end

% prepare general Parameter structure P
P.file_date  = DateTime;
P.file_dur = Duration;
P.Fs_SP = Freq;
P.Fs_AD = Fs;
P.Gains_SP = gains_ap;
P.Gains_AD = gains_ad;
P.AD_enabled = ad_chan_on;

% ----------- first convert data with stimulation -------------------------            
if spontaneous == 0

    % convert header in MATLAB-struct array
    code_start = find(sv == 253,1,'first'); code_stop = find(sv == 253,1,'last');
    % if there is only one 253 test whether data are intact and add one 253
    if code_start == code_stop
        if all(sv(1:code_stop-1) < 253 & (sv(1:code_stop-1) >= 230 | sv(1:code_stop-1) == 0))
         sv = [253;sv];
         evts = [0;evts];
         nev = [nev+1];
         code_start = 1; code_stop = code_stop+1;
        end
    end
    code       = sv(code_start:code_stop);
    % take out zeros these are AT LEAST IN THE HEADER not used at all
    code = code(find(code > 0));
    
    % due to bad behavior of events take out zeros during the code
    try
        [Header] = Online_STRBCode2Header_01(code);
        % new code to access m-files
        header_fields = fieldnames(Header);
        if all(isequal(header_fields{:},'file'))
            p_cur = cd;
            hostname = '\\pcr5510-4';
            n_exp_files = [hostname,Header.file(1,strfind(Header.file,':')+1:end)];            
            [p_exp_file,n_exp_file] = fileparts(n_exp_files);
            p_parent_file = [fileparts(fileparts(p)),filesep];
            % take the associated m-file from plx folder if it exists
            if exist([p_parent_file,n_exp_file,'.m'],'file')
                cd(p_parent_file);
                n_exp_files = [p,n_exp_file,'.m'];
            elseif exist(Header.file,'file')
                n_exp_files = [Header.file(1,strfind(Header.file,':')+1:end)];
                [p_exp_file,n_exp_file] = fileparts(n_exp_files);
                cd(p_exp_file)
            elseif exist([p,n_exp_file,'.m'],'file')
                cd(p);
                n_exp_files = [p,n_exp_file,'.m'];
            else
                cd(p_exp_file);  
            end
            Header = feval(n_exp_file);
            Header.file = n_exp_files;
            cd(p_cur);
        end
        Header.ncond = size(Header.stimlist,1);
        
        % take out zeros these are WRONG strobes 0 is not used at all
        i_sv_ok = find(sv > 0);
        sv = sv(i_sv_ok);
        evts = evts(i_sv_ok);
    catch
        Header = Play_Code2Header_01(code); % Matthias setup
        try
        Header = Header(1);
%         catch
            keyboard;
        end
        % convert some variables into a common scheme
        Header.n_rep = Header.nrep;
        Header.ISI   = [Header.iti_fix Header.iti_var+Header.iti_fix];
        Header.dura  = 0; % just a dummy for now
        un_stim = unique(Header.s2);
        dummy_rep = zeros(size(Header.s2,1),1);
        for i1 = 1:size(un_stim,1)
            i_stim = find(Header.s2 == un_stim(i1));
            dummy_rep(i_stim) = [1:size(i_stim,1)];
        end
        Header.evtlist  = [Header.s2,dummy_rep];
        Header.stimlist = cat(2,Header.s1',mat2cell(zeros(size(Header.s1',1),2),ones(size(Header.s1',1),1),ones(2,1)));
        Header.exp_type = 'Matthias_311';
        % for Matthias setup do not take out strobes 0
    end
    if isempty(Header)
        %display notification about found strobe events
        disp(sprintf('Header is missing file will not be converted!'));
        return;
    end
    Header = Header(1);

    % recheck input here because in Header stimulus durations should be
    % given
    if dur - t_pre <= max(Header.dura)
        dur = max(Header.dura)+ t_pre*2;
    end
    if dur - t_pre >= min(Header.ISI)
        dur = min(Header.ISI)+t_pre;
    end
    if max(Header.dura) > min(Header.ISI) & max(Header.dura)+t_pre <= max(Header.dura) + min(Header.ISI)
        dur = max(Header.dura)+t_pre*2;
    elseif max(Header.dura) > min(Header.ISI) & max(Header.dura)+t_pre > max(Header.dura) + min(Header.ISI)
        dur = min(Header.ISI)+max(Header.dura)+t_pre;
    end

    
    stim_names      = unique(sv(find(sv < 230))); n_stim = length(stim_names); n_pres = length(sv(find(sv < 230)));% line not finished yet what is upper limit for number of stimuli (230?)?
    stim_list       = sv(find(sv(code_stop+1:end) < 230)+code_stop);
    % list of timestamps of stimuli
    tts = evts(find(sv(code_stop+1:end) < 230)+code_stop);

    % check whether number of events are the same as expected
    if size(Header.n_rep,2) == Header.ncond
        n_exp_evts = sum(Header.n_rep);
    else
        n_exp_evts = Header.ncond*Header.n_rep;
    end
    if ~(size(tts,1) == n_exp_evts)
        % if not check whether number of trigger are same as number of
        % strobes cause this would indicate that user stopped recording
        % earlier but file integrity is still given
        if ~(n_trig == size(tts,1))
             %display notification about found strobe events
             disp(sprintf('Inconsistent number of trigger %d and strobes %d!',n_trig, size(tts,1)));
             disp(sprintf('YOUR FILE IS MOST LIKELY SCREWED.'));
             disp(sprintf('Attempting work around...'));
                 
            % all right here new code has to be added to deal with that   
            % check whether some strobes followed to close to each other
            tts_diff = diff(tts);
            i_spurious = find(tts_diff <= 0.01); %strobes faster then 10ms
            
            % now check whether those spurious items are the same strobe
            ind_notok = [];
            for i1 = 1:size(i_spurious,1)
                if stim_list(i_spurious(i1)) == stim_list(i_spurious(i1)+1)
                    ind_notok(1,end+1) = i_spurious(i1)+1;
                end
            end
            ind_ok = setxor(ind_notok,[1:size(tts,1)]);
            
            trial_ts = tts;
             
            %display notification about found strobe events
            disp(sprintf('Found %d strobes that are too close and\nhave the same identity.',size(ind_notok,2)));
            disp(sprintf('That makes %d strobes and %d expected events!',size(ind_ok,2),n_exp_evts));
            if n_exp_evts == size(ind_ok,2)
                disp(sprintf('YOUR PROBLEM HAS BEEN SOLVED.'));
            else
                disp(sprintf('YOUR PROBLEM COULD NOT BE SOLVED.'));
            end
        else
            trial_ts = ts_trig; % change from strobe timing to trigger timing
            ind_ok = [1:size(tts,1)];
        end
    else 
        if n_trig == size(tts,1)
            trial_ts = ts_trig; % change from strobe timing to trigger timing
        else
            trial_ts = tts;
        end
        ind_ok = [1:size(tts,1)];
    end 
    % correct stimulus list and list of timestamps of stimuli
    stim_list = stim_list(ind_ok);
    trial_ts = trial_ts(ind_ok);
    n_pres = size(stim_list,1);     

    % allow to select certain events
    EvtIdentifier = '';
    % THIS IS FOR CHOOSING SUBSET OF STIMULI uncommented by Max 17.06.2015
    %answer = inputdlg(sprintf('Found %i events for %i stimuli. Which ones to take?',n_pres,Header.ncond),'Event select',1,{sprintf('1, %i',n_pres)});
    answer=[1 n_pres];
    try
        evt_select = str2num(answer{1});   
        if length(evt_select) == 1
            if evt_select(1) < 0
                evt_select = [length(stim_list)+evt_select+1,length(stim_list)];
            else evt_select(1) > 0
                evt_select = [1, evt_select];
            end
        end
        if evt_select(2)-evt_select(1)+1 ~= size(stim_list,1);
            EvtIdentifier = sprintf('_Evt%i-%i',evt_select);
        end
        stim_list = stim_list(evt_select(1):evt_select(2));
        trial_ts = trial_ts(evt_select(1):evt_select(2));  
        n_pres   = size(stim_list,1);
    end
    
    % Read AD data
    if ~isempty(n_ad_chan_on) & n_ad_chan_on ~= 0;
        % Create t vector
        t_AD = [1/Fs:1/Fs:dur];
        n_ADpoints = length(t_AD);
        V=zeros(n_ad_chan_on,n_ADpoints,n_pres);
        
        hdl_WaitBar = waitbar(0, 'Reading AD data ...');
        c_sweep = 1;

        for c = 1:n_ad_chan_on
            [adfreq, n_pointstotal, ts, fn, AD] = mexPlex(17,[p,n], ad_chan_on(c)-1);
            if length(ts) ~= 1
                        disp(sprintf('WARNING : AD data (channel index %d) fully reconstructed from %d fragments.', 1, length(fn)));        
            end
            for r = 1 : n_pres   
                    i0 = ceil((trial_ts(r) - ts(1)-t_pre) * Fs);
                    i1 = i0 + n_ADpoints -1;
                    if length(AD) >= i1
                        V(c,:,r) = AD(i0:i0+n_ADpoints-1);
                    else
                        disp(sprintf('WARNING : AD data (chan.ind.=%d, r=%d,s=%d) incomplete.', c, r, trial_ts(r)));        
                        V(c,:,r) = zeros(n_ADpoints,1);
                    end  
                waitbar(c_sweep / (n_ad_chan_on*n_pres));
                c_sweep = c_sweep+1;
            end
        end     
        close(hdl_WaitBar);
        drawnow;

        AD = V;
    else
        disp('No AD Data found!');
        AD = [];
    end

    % Read AP data
    AP  = [];   % dimensions of Spikes-Matrix: r: occurence; c: 1 - trial, 
                        % 2 - stimulus, 3 - ID-channel, 4 - ID-unit, 5 - ts 
                        % according to file start, 6 - ts according to stimulus-
                        % onset, 7-38 - waveform of this spike (if all zeros - no
                        % waveform recorded
       
    if ~isempty(ID)
        disp('Spike channel and unit IDs :');
        disp(ID');
        hdl_WaitBar = waitbar(0, 'Reading Unit data ...');
        drawnow;
        c_sweep = 1;

        for c = 1:n_un_tot
             try 
%                 [n_spik_un, npw, spik_ts, wave] = mexPlex(19,[p,n], ID(c,1), ID(c,2)-1);
%             catch
                [n_spik_un, spik_ts] = mexPlex(5, [p,n], ID(c,1), ID(c,2)-1);
                wave = zeros(n_spik_un,NPW);
            end
            for r = 1 : n_pres   
                    t1 = trial_ts(r)-t_pre;
                    t2 = t1 + dur;
                    spik_ind     = find(spik_ts(:,1) >= t1 & spik_ts(:,1) <= t2);
                    n_spik_ind   = length(spik_ind);
                    dummy_spik   = [];
                    if ~isempty(n_spik_ind)
                        dummy_spik   = [repmat([r,stim_list(r),ID(c,1),ID(c,2)],[n_spik_ind,1]),spik_ts(spik_ind,1),...
                            spik_ts(spik_ind,:)-trial_ts(r),wave(spik_ind,:)];
                        AP = [AP;dummy_spik];
                    end      
                    waitbar(c_sweep / (n_un_tot*n_pres));
                c_sweep = c_sweep+1;
            end            
        end     
        close(hdl_WaitBar);
        drawnow;
    else
        disp('No Spike Channel found!');
        AP = [];
    end

    % Create Spik_list and avgFP for OnlineEval-tools
    [Spik_list, avgFP, n_avg] = Create_Spiklist_AvgFP_00(Header,stim_list,AP,AD);

    % update the Header 
%     if iscell(Header.stimlist)
%         Header.t_sig = [Header.stimlist{:,end-1}];
%     else
%         Header.t_sig = Header.stimlist(:,end-1)';
%     end
    Header.t_sig = Header.dura';
    Header.t_pre = t_pre;  Header.t_pos = dur-Header.t_sig-t_pre;
    
    % figure out how many reps per stim are actually in this file
    un_stim   = unique(stim_list);
    un_stim   = un_stim(find(~isnan(un_stim)));
    n_stim    = length(un_stim);
    stim_ind  = zeros(size(un_stim));

    n_stim_rep = 0;
    for i1 = 1:size(un_stim,1)
         stim_ind(i1) = size(find(stim_list == un_stim(i1)),1);
    end
    Header.n_rep = stim_ind';
    if all(Header.n_rep(1) == Header.n_rep)
        Header.n_rep = Header.n_rep(1);
    end

    % save data
    if type_choi == 1   & save_choi == 1
            DS  = Convert_Header2SWEEP_01(Header,t_pre,dur,stim_list);

            % put everything into the SWEEP-structure
            % Init SWEEP structure
            SWEEP(size(stim_list,1),1) = struct('Header',[], 't0', [], 'Spikes',[], 'AD',[]);

            c_event2 = 1;
            for r = 1 : size(stim_list,1)
                        SWEEP(r).Header = stim_list(r);
                        SWEEP(r).t0     = trial_ts(r);
                        if ~isempty(AD)
                            SWEEP(r).AD     = AD(:,:,r)';
                        end
                        if ~isempty(AP)
                        SWEEP(r).Spikes = AP(find(AP(:,1) == r),[6 3 4 7:end]);
                        end
                    c_event2 = c_event2+1;
            end
            if size(DS.n_Rep,2) == 1 & size(DS.n_Stm,2) == 1 & size(SWEEP,1) == size(DS.n_Stm,1)*size(DS.n_Rep,1)
                SWEEP           = reshape(SWEEP,[DS.n_Stm,DS.n_Rep])';
            end
            
            disp(        'Saving!');
            save([p,n(1,1:end-4),EvtIdentifier,'.mat'], 'DS', 'P', 'SWEEP', 'ID', 'Header','Spik_list','avgFP'); % modified by Max / 17.06.2015
            save_choi = 0;
    end  
    
% ---------------- now convert recordings of spontaneous activity ---------
else spontaneous == 1 
    
    disp('Spontaneous file encountered!');
    [DS,Header,AD,AP] = Convert_SpontPLX2MAT_00(p,n,ID,ad_chan_on);

    % create last missing output variable/s
    Header.animal_ID = n(1,1:min(strfind(n,'_'))-1);
    Header.datetime  = datestr(datenum(DateTime,'mm/dd/yyyy HH:MM:SS'),'yyyy-mm-dd HH:MM:SS');
    
    if type_choi == 1  & save_choi == 1
            SWEEP(1).AD = AD;
            SWEEP(1).Spikes = AP;
            SWEEP(1).Header = [];
            
            disp(        'Saving!');
            save([p,n(1,1:end-4),'.mat'], 'DS', 'P', 'SWEEP', 'ID', 'Header');%
            save_choi = 0;
    end
end

% ---------------- save if user has chosen to -----------------------------
if save_choi == 1
    disp(        'Saving!');
    if ~isempty(AD)
        save([p,n(1,1:end-4),'_AD.mat'],'AD');
    end
    if ~isempty(AP)
        save([p,n(1,1:end-4),'_AP.mat'],'AP');
    end 
end

% COMMENTED ON 2012-07-12 BECAUSE FILE IS NOT NECESSARY
% % ---------------- always save parameter file -----------------------------
% if type_choi == 1  | save_choi == 1
%     save([p,n(1,1:end-4),'_PAR.mat'],'DS','ID','P','Header');
% end
disp(        'Finished!');
disp(        '*************************************************');
clear all
