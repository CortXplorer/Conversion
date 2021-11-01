function  [DATA,P,ID] = Convert_NISetup_PLX2MAT_simple_10(filename,save_choi,type_choi,t_pre,dur, TLE,flag_spk);
% Should be used like this: 
% prestimulusinterval = 0;
% poststimulusinterval = duration of recording session or inf;
% durationofepoch =prestimulusinterval+poststimulusinterval;
% flag_spk = 1 if spike output required;
% 
% [dat,P,ID] = Convert_NISetup_PLX2MAT_simple_10(pathfilename,0,0, prestimulusinterval, durationofepoch ,0,flag_spk);



% Convertion here safes values in V!
% check input

if ~exist('filename') || isempty(filename)
    pth_data = 'D:\MyCode\Conversion';
    pth_src = pwd;
    cd(pth_data);
    filename = '';
    [n,p] = uigetfile('*.plx');
    cd(pth_src);
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
if ~exist('TLE') | isempty(TLE)
    TLE   = [NaN];
end

events0 = struct('label',[],'chan',[],'code',[],'datenum',[],'t',[],'ts',[],'dur',[],'order',[],'ts_epoch',[],'order_epoch',[],'iepoch',[],'session_file',[]);
epochs0 = struct('label',[],'chan',[],'code',[],'Ts',[],'T0',[],'datenum',[],'t',[],'ts',[],'ievent',[],'iepoch',[],'session_file',[],'events',[]);


% stop execution if no file has been chosen
if n == 0; return; end

% inform about file name
disp(sprintf('\n'));
disp('*************************************************');
disp(sprintf('PROCESSING          : %s',n));
disp(sprintf('IN FOLDER           : %s\n',p));


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
try
    [nev, evts, sv] = mexPlex(3,[p,n], 257);
catch
    nev = 0; evts= 0; sv = 0;
end

events = repmat(events0,[1,nev+n_trig]);      
for ie =1:nev
    events(ie) = struct('label',[num2str(sv(ie))],'chan',[],'code',[sv(ie)],'datenum',[],'t',[evts(ie)],'ts',[evts(ie)],'dur',[0],'order',[NaN],'ts_epoch',[NaN],'order_epoch',[NaN],'iepoch',[NaN],'session_file',[1]);
end
for ie =1:n_trig
    events(nev+ie) = struct('label',[num2str(257)],'chan',[0],'code',[257],'datenum',[],'t',[ts_trig(ie)],'ts',[ts_trig(ie)],'dur',[0],'order',[NaN],'ts_epoch',[NaN],'order_epoch',[NaN],'iepoch',[NaN],'session_file',[1]);
end
ts = [events(:).ts];
[j,ind] = sort(ts(:));
events = events(ind);

% prepare general Parameter structure P
P.file_date  = DateTime;
P.file_dur = Duration;
P.Fs_SP = Freq;
P.Fs_AD = Fs;
P.Gains_SP = gains_ap;
P.Gains_AD = gains_ad;
P.AD_enabled = ad_chan_on;
% ----------- EPOCHS -------------------------            
if isnan(TLE)
   TLE = unique(sv(:));
end
ind = [];
for itle = 1:length(TLE(:))
    iii = find(sv(:) == TLE(itle));
    ind = [ind(:);iii(:)];
end

code = sv(ind);
tsepoch = evts(ind);
[tsepoch,ind] = sort(tsepoch(:));
code = code(ind);
nepoch = length(tsepoch(:));


epochs = repmat(epochs0,[1,nepoch]);
for ie = 1:nepoch
    iev = find([events(:).ts] == tsepoch(ie));
    ievent = find([events(:).ts] >= tsepoch(ie)-t_pre & [events(:).ts] < tsepoch(ie)+dur-t_pre);
    epochs(ie) = struct('label',[events(iev).label],'chan',[events(iev).chan],'code',[events(iev).code],'Ts',[-t_pre + dur-t_pre],'T0',[-t_pre + dur-t_pre],...
                        'datenum',[events(iev).datenum],'t',[events(iev).t],'ts',[events(iev).ts],'ievent',[iev],'iepoch',[ie],'session_file',[NaN],'events',[events(ievent)]);
    [events(ievent).iepoch] = deal(ie);               
end

% ----------- first convert data with stimulation -------------------------            
% Read AD data
if ~isempty(n_ad_chan_on) && n_ad_chan_on ~= 0;
    % Create t vector
    t_AD = [1/Fs:1/Fs:dur];
    n_ADpoints = length(t_AD);
    V=zeros(n_ad_chan_on,n_ADpoints,nepoch);
    hdl_WaitBar = waitbar(0, 'Reading AD data ...');
    c_sweep = 1;

    for c = 1:n_ad_chan_on
        [adfreq, n_pointstotal, ts, fn, AD] = mexPlex(17,[p,n], ad_chan_on(c)-1);
        if length(ts) ~= 1
                    disp(sprintf('WARNING : AD data (channel index %d) fully reconstructed from %d fragments.', 1, length(fn)));        
        end
        for r = 1 : nepoch  
                i0 = ceil((tsepoch(r) - ts(1)-t_pre) * Fs);
                i1 = i0 + n_ADpoints -1;
                if length(AD) >= i1
                    V(c,:,r) = AD(i0:i0+n_ADpoints-1);
                else
                    disp(sprintf('WARNING : AD data (chan.ind.=%d, r=%d,s=%d) incomplete.', c, r, tsepoch(r)));        
                    V(c,:,r) = zeros(n_ADpoints,1);
                end  
            waitbar(c_sweep / (n_ad_chan_on*nepoch));
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
       
if ~isempty(ID) && flag_spk == 1
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
        for r = 1 : nepoch   
                t1 = tsepoch(r)-t_pre;
                t2 = t1 + dur;
                spik_ind     = find(spik_ts(:,1) >= t1 & spik_ts(:,1) <= t2);
                n_spik_ind   = length(spik_ind);
                dummy_spik   = [];
                if ~isempty(n_spik_ind)
                    dummy_spik   = [repmat([r,code(r),ID(c,1),ID(c,2)],[n_spik_ind,1]),spik_ts(spik_ind,1),...
                        spik_ts(spik_ind,:)-tsepoch(r),wave(spik_ind,:)];
                    AP = [AP;dummy_spik];
                end      
                waitbar(c_sweep / (n_un_tot*nepoch));
            c_sweep = c_sweep+1;
        end            
    end     
    close(hdl_WaitBar);
    drawnow;
else
    disp('No Spike Channel found!');
    AP = [];
end



%cd(pth_data);
   % save data
if type_choi == 1   & save_choi >= 1

        % put everything into the SWEEP-structure
        % Init SWEEP structure
        SWEEP(size(code,1),1) = struct('Header',[], 't0', [], 'Spikes',[], 'AD',[],'code',[]);

        c_event2 = 1;
        for r = 1 : nepoch
                    SWEEP(r).Header = code(r);
                    SWEEP(r).t0     = tsepoch(r);
                    SWEEP(r).code     = code(r);
                    if ~isempty(AD)
                        SWEEP(r).AD     = AD(:,:,r)';
                    end
                    if ~isempty(AP)
                    SWEEP(r).Spikes = AP(find(AP(:,1) == r),[6 3 4 7:end]);
                    end
                c_event2 = c_event2+1;
        end

        disp(        'Saving!');
       save([p,n(1,1:end-4),'.mat'], 'P', 'SWEEP', 'ID'); % modified by Max / 17.06.2015
      
else
     
   DATA.LFP =  AD; %channels x time x epoch 
   DATA.SPK = [];
   if ~isempty(AP)
       DATA.SPK =AP(find(AP(:,1) == r),[6 3 4 7:end]);
   end
   DATA.EPOCHS = epochs;
   DATA.EVENTS = [events];
   DATA.ART = [];
   DATA.TYPE = [];
   DATA.PATH = '';
   if  ~isempty(DATA.PATH) && ~isequal(DATA.PATH(end),'\')
       DATA.PATH(end) = '\';
   end
   DATA.FILE = n;
   DATA.FILE_VID = [];
   
   
   DATA.DATE = DateTime;
   DATA.SR = P.Fs_AD(1);
   DATA.TLIM =  [-t_pre -t_pre + dur];
   DATA.CH = [];
   DATA.CHMAP = [];
   DATA.CHNAMES = [];
   DATA.ELECPOS = [];

   AD = [];
   AP = [];
   
end  

% ---------------- save if user has chosen to -----------------------------

if save_choi >= 2
    disp(        'Saving!');
    
    if ~isempty(AD)
        save([p,n(1,1:end-4),'_AD.mat'],'AD');
    end
    if ~isempty(AP)
        save([p,n(1,1:end-4),'_AP.mat'],'AP');
    end 
    if isempty(AD) & isempty(AP)
         save([p,n(1,1:end-4),'_AP.mat'],'DATA');
    end
end



% COMMENTED ON 2012-07-12 BECAUSE FILE IS NOT NECESSARY
% % ---------------- always save parameter file -----------------------------
% if type_choi == 1  | save_choi == 1
%     save([p,n(1,1:end-4),'_PAR.mat'],'DS','ID','P','Header');
% end
disp(        'Finished!');
disp(        '*************************************************');
%clear all
%   cd(pth_src);