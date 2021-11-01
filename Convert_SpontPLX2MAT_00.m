function [DS,Header,AD,AP] = Convert_SpontPLX2MAT_00(p,n,ID,ad_chan_on)

%% check input
if ~exist('p') | strmatch(p,'') | ~exist('n') | strmatch(n,'')
    [n,p] = uigetfile('*.plx');
    % Calc file size of PLX file
    f = dir(p);
    size_file = f(strmatch(n,char(f.name))).bytes;
    disp(sprintf('SIZE OF FILE          : %.1f KB, %.2f MB', size_file/1024, size_file/1024/1024));
    
     % read general data from the plx-file
    % general
    [OpenedFileName, Version, Freq, Comment, Trodalness, NPW, PreTresh, SpikePeakV, SpikeADResBits, SlowPeakV, SlowADResBits, Duration, DateTime] = mexPlex(13,[p,n]);
    
    % AP-Data
    [n_chan_tot,thresh]  = mexPlex(9,[p,n]);   % output is array of threshold vals for each channel
    [n_chan_tot,gains_ap]   = mexPlex(8,[p,n]);    % output is total gain (including preamp gain and digital gain)
    [n_chan_tot,names]   = mexPlex(14,[p,n]);
    [n_chan_tot,filters] = mexPlex(10,[p,n]); % check which spikechannels have been filtered
    % preliminary way to check which channels had been enabled (write email
    % to plexonincs
    mingain = 1000;
    ap_chan_on = find(gains_ap > mingain);n_ap_chan_on = length(ap_chan_on);
    % get information about number of units, number of spikes ...
    [tscounts, wfcounts, evcounts] = mexPlex(4,[p,n], 1); % get info about enabled channels and units

    [dummy_row, dummy_col] = find(tscounts);
    ID = [dummy_col-1, dummy_row-1];
    n_un_tot = size(ID,1);
    
    % AD-Data
    [n_chan_tot,freqs] = mexPlex(12,[p,n]); % output is sampling frequency; n_chan is number of channels the system supports
    [n_chan_tot,gains_ad] = mexPlex(11,[p,n]); % output is total gain (including preamp gain and digital gain)
    [n_chan_tot,names] = mexPlex(15,[p,n]);
    % preliminary way to check which channels had been enabled (write email
    % to plexoninc
    preampgain = 1000;
    ad_chan_on = find(gains_ad > preampgain);n_ad_chan_on = length(ad_chan_on);
    Fs =freqs(1);  
    
    % prepare general Parameter structure P
    P.file_date  = DateTime;
    P.file_dur = Duration;
    P.Fs_SP = Freq;
    P.Fs_AD = Fs;
    P.Gains_SP = gains_ap;
    P.Gains_AD = gains_ad;
    P.AD_enabled = ad_chan_on;
end

%% get necessary data for SPikes
if ~exist('ID') | isempty(ID)
    % AP-Data
    [n_chan_tot,thresh]  = mexPlex(9,[p,n]);   % output is array of threshold vals for each channel
    [n_chan_tot,gains_ap]   = mexPlex(8,[p,n]);    % output is total gain (including preamp gain and digital gain)
    [n_chan_tot,names]   = mexPlex(14,[p,n]);
    [n_chan_tot,filters] = mexPlex(10,[p,n]); % check which spikechannels have been filtered
    % preliminary way to check which channels had been enabled (write email
    % to plexonincs
    mingain = 1000;
    ap_chan_on = find(gains_ap > mingain);n_ap_chan_on = length(ap_chan_on);
    % get information about number of units, number of spikes ...
    [tscounts, wfcounts, evcounts] = mexPlex(4,[p,n], 1); % get info about enabled channels and units

    [dummy_row, dummy_col] = find(tscounts);
    ID = [dummy_col-1, dummy_row-1];
end
n_un_tot = size(ID,1);

%% get necessary data for Field Potentials
if ~exist('ad_chan_on') | isempty(ad_chan_on)
    % AD-Data
    [n_chan_tot,freqs] = mexPlex(12,[p,n]); % output is sampling frequency; n_chan is number of channels the system supports
    [n_chan_tot,gains_ad] = mexPlex(11,[p,n]); % output is total gain (including preamp gain and digital gain)
    [n_chan_tot,names] = mexPlex(15,[p,n]);
    % preliminary way to check which channels had been enabled (write email
    % to plexoninc
    preampgain = 1000;
    ad_chan_on = find(gains_ad > preampgain);
    Fs =freqs(1);
end
n_ad_chan_on = length(ad_chan_on);

%% Read AP data
if ~isempty(ID)
    disp('Spike channel and unit IDs :');
    disp(ID');
    hdl_WaitBar = waitbar(0, 'Reading Unit data ...');
    drawnow;
    c_sweep = 1;AP = [];
    for c = 1:n_un_tot
        try
            [n_spik_un, npw, spik_ts, wave] = mexPlex(19,[p,n], ID(c,1), ID(c,2));
        catch
            [n_spik_un, spik_ts] = mexPlex(5, [p,n], ID(c,1), ID(c,2));
            wave = zeros(n_spik_un,NPW);
        end
        dummy_spik   = [repmat([1,1,ID(c,1),ID(c,2)],[size(spik_ts,1),1]),spik_ts,...
            spik_ts,wave];
        
        % added this fix for sizing/concatonation issues! KD - 08-10-2019
        if c > 1 && size(dummy_spik,2) > size(AP,2)
            AP = horzcat(AP,nan(size(AP,1),size(dummy_spik,2)-size(AP,2)));
        end
        if c > 1 && size(dummy_spik,2) < size(AP,2)
            dummy_spik = horzcat(dummy_spik, nan(1,size(AP,2)-size(dummy_spik,2)));
        end
        
        AP = [AP;dummy_spik];
        waitbar(c_sweep / (n_un_tot));
        c_sweep = c_sweep+1;
    end
    close(hdl_WaitBar);
    drawnow;
else
    disp('No Spike Channel found!');
    AP = [];
end

%% Read AD data
if ~isempty(n_ad_chan_on)
    hdl_WaitBar = waitbar(0, 'Reading AD data ...');
    c_sweep = 1;

    for c = 1:n_ad_chan_on
        [adfreq, n_pointstotal, ts, fn, AD] = mexPlex(17,[p,n], ad_chan_on(c)-1);
        if length(ts) ~= 1
                    disp(sprintf('WARNING : AD data (channel index %d) fully reconstructed from %d fragments.', 1, length(fn)));        
        end
        V(c,:,1) = AD;

        waitbar(c_sweep / (n_ad_chan_on));
        c_sweep = c_sweep+1;
    end     
    close(hdl_WaitBar);
    drawnow;

    AD = V;
else
    disp('No AD Data found!');
    AD = [];
end
    
%% create some simple entrys in DS and Header structure
DS.type          = 0;
DS.n_Rep		 = 1;
DS.n_Stm         = 1;
DS.t_Pre         = 0;
DS.t_Sig         = 0;
DS.t_Pau         = size(AD,2)/adfreq; % takes length of spontaneous recording as time

Header.type      = 'spontaneous';
Header.code_type = 0;
Header.nrep      = 1;
Header.ncond     = 1;
Header.t_pre     = 0;
Header.t_sig     = 0;
Header.t_pos     = size(AD,2)/adfreq; % takes length of spontaneous recording as time
