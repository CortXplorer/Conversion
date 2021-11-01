function varargout = Play_ParaCheck_02(varargin)

%input
Hw.nio1.chans = [0:7]; 

TypePara = {};
PP = {}; HH = {}; FF = {};
ntype = 0; npara = 0;
count_type = 0; count_in = 0;
while count_in < nargin
    count_in = count_in + 1;
    if isstr(varargin{count_in})
       count_type = count_type + 1;
       TypePara{count_type} = varargin{count_in};
       PP{count_type} = [];
       HH{count_type} = [];
    end
    if isstruct(varargin{count_in})
       PP{count_type} = varargin{count_in};
       if isequal(lower(TypePara{count_type}), 'play'); 
          count_in = count_in + 1;
          HH{count_type} = varargin{count_in}; 
       else; 
          HH{count_type} = []; 
       end;
    end;
    if isnumeric(varargin{count_in}) & count_in == nargin
       Hw = varargin{nargin};
    end
end

%Defaults
for it = 1:count_type
    if isempty(PP{it}); flag_def = 0; else; flag_def = 1; end;
    switch lower(TypePara{it})
        case 'play'

            H0   = struct('type',['play'],'type_play',['none'],'code_play',[],'info',[''],'header',[],'code_header',[],'list_play',[''],'lpt1',[],'nio1',[],'nid1',[],'stg1',[],'com1',[],...
                             'Ind',[],'SubInd',[],'Isi',[],'isi_fix',[0],'isi_var',[inf]);            
            H0.lpt1  = struct('type',['lpt1'],'status',0,'name',['parallelLPT1-DIO'],'device',['1'],'id',[''],'errmsg',[]); 
            H0.nio1  = struct('type',['nio1'],'status',0,'name',['nidaq1-AO'],'device',[''],'id',['1'],'nch',[''],'errmsg',[]);
            H0.nid1  = struct('type',['nid1'],'status',0,'name',['nidaq1-DIO '],'device',[''],'id',['1'],'mode','in','code',[],'errmsg',[]); 
            H0.stg1  = struct('type',['stg1'],'status',0,'name',['STG-Isolator-8'],'device',[''],'id',['1'],'nch',[8, 4],'errmsg',[]);
            H0.com1  = struct('type',['com1'],'status',0,'name',['g.PAH'],'device',[''],'id',['1'],'code0',[4*16+13 5*16+10 6*16+9],'errmsg',[]);
            F0_H  = [fieldnames(H0),  {'str';'str';'num';'str';'stc';'numarray';'cellstr';'sub';'sub';'sub';'sub';'sub';'numarray';'numarray';'numarray';'num';'num'}];
            F1_H  = {'nio1';'nid1';'lpt1';'stg1';'com1'};

            
            P0       = struct('type',['play'],'nio1',[],'nid1',[],'lpt1',[],'stg1',[],'com1',[],'info',['']);
            P0.nio1  = struct('type',['nio1'],'data',[],'chan',[],'range',[-10 10],'nrep',[0],'type_trig', ['HwDigital'], 'flag_off', [0], 'sr',[],'nch',0,'npt',0,'status',[0]);
            P0.nid1  = struct('type',['nid1'],'status',[0]);
            P0.stg1  = struct('type',['stg1'],'data',[],'t',[],'chan',[],'mode', ['I'], 'nrep',[0],'type_trig', [''], 'sr',[10^6],'nch',0,'npt',[0],'status',[0]);
            P0.lpt1  = struct('type',['lpt1'],'code',[0 0 0 0],'status',[0]);
            P0.com1  = struct('type',['com1'],'code',[],'status',[0]);
            F0_P  = [fieldnames(P0),  {'str';'sub';'sub';'sub';'sub';'sub';'str'}];
            F1_P  = {'nio1';'nid1';'lpt1';'stg1';'com1'};

        case 'tone1'
            H0 = [];
            F0_H  = {};
            F1_H  = {};
            callib = [0 0; 250 0; 500 0; 1000 0; 2000 0; 4000 0; 8000 0; 16000 0; -inf 0; inf 0]';
            P0 = struct('type', ['tone1'], 'code_type',[101], 'f',[],'att',[0],'dur',[200],'ph',[0],'dur_click',[0.1],'dur_ramp',[5],'type_ramp',['linear'],'sr',[120000],'delay_pre',[0],'srep',[0],'nrep',[100],...
                        'isi_fix',[1.0],'isi_var',[0],'order',[0 0 0],'rand',[1 0 0], 'callib',callib,'chan_sig', Hw.nio1.chans(1), 'range_sig', [-3 3], 'type_trig','HwDigital' ,'chan_trig', Hw.nio1.chans([3,5]),'ncond',0,'Cond',[],'var_cond',[],'COND',[]);
            F0_P  = [fieldnames(P0),  {'str';'num';'numarray';'numarray';'numarray';'num';'num';'num';'str';'num';'num';'num';'num';'num';'num';'numarray';'numarray';'numarray';'numarray';'numarray'; 'str';'numarray';'num';'struct';'cellarray';'cellarray'}];
            F1_P  = {};

        case 'puls1'
            H0 = [];
            F0_H  = {};
            F1_H  = {};
            P0 = struct('type', ['puls1'], 'code_type',[901],'a_ph1',[],'chan',[],'dur_ph1',[200],'range',[-800, 800], 'iphi',[100],'scale_ph2',[-1],'gate_puls',[-1.5 1.5],'sr',[50000],'delay_pre',[0],'srep',[0],'nrep',[100],...
                        'isi_fix',[1.0],'isi_var',[0],'order',[0],'rand',[1 0 0], 'type_trig','HwDigital' ,'chan_trig', Hw.nio1.chans([4,5]),'chan_sync', [0,3], 'ncond',0,'Cond',[],'var_cond',[],'COND',[]);
            F0_P  = [fieldnames(P0),  {'str';'num';'numarray';'numarray';'numarray';'numarray'; 'num';'num';'numarray';'num';'num';'num';'num';'num';'num';'numarray';'numarray';'str'; 'numarray';'numarray';'num';'struct';'cellarray';'cellarray'}];
            F1_P  = {};
            
        case 'cftse'
            H0 = [];
            F0_H  = {};
            F1_H  = {};
            P0 = struct('type','cftse','code_type',[001],'file_ctr','','type_ctr','','iti_fix',[],'iti_var',[],'wait_trial',[],'us_min',[],'us_max',[],'liampl',[],'hurdle',[],'hurdledelay',[],'hurdleampl',[],'s1',[],'s2',[],'nrep',[],'ncond',[],'Cond',[],'var_cond',[],'COND',[]);
            F0_P = [fieldnames(P0), {'str';'num';'str';'str';'num';'num';'num';'num';'num';'num';'str';'num';'num';'cellstr';'numarray';'num';'num';'struct';'cellarray';'cellarray'}];
            F1_P  = {};

        case 'fmsweep'
            H0 = [];
            F0_H  = {};
            F1_H  = {};
            P0 = struct('type', ['fmsweep'], 'code_type',[102], 'f_start',[0],'range_f',[10000],'f_mod',[2],'type_mod',['linear full cycle'],'att',[20],'dur_ramp',[5],'type_ramp',['linear'],'sr',[120000],'nrep',[100],...
                        'chan_sig', Hw.nio1.chans(1), 'range_sig', [-3 3], 'type_trig','HwDigital' ,'chan_trig', Hw.nio1.chans([3,4]),'ncond',0,'Cond',[],'var_cond',[],'COND',[]);
            F0_P  = [fieldnames(P0),  {'str';'num';'num';'num';'num';'str';'num';'num';'str';'num';'num';'numarray';'numarray';'str';'numarray';'num';'struct';'cellarray';'cellarray'}];
            F1_P  = {};
            
        case 'FM'
        case 'AM'
            
        otherwise
            H0 = []; P0 = []; F0 = {}; F1 = {};
    end

    %Field Check
    if flag_def == 1; 
       for ihp = 1:2 
           if ihp == 1; S = PP{it};  S0 = P0; F0 = F0_P; F1 = F1_P; end;
           if ihp == 2; S = HH{it};  S0 = H0; F0 = F0_H; F1 = F1_H; end;
           n_S = length(S(:));
           
           for i = 1:n_S
                %check field existence and set default values
                for ifn = 1:size(F0,1)
                    if ~isfield(S(i), F0{ifn,1}) | isempty(S(i).(F0{ifn,1}))
                       S(i).(F0{ifn,1}) = S0.(F0{ifn,1});
                    end
                    for ifn_sub = 1:length(F1(:))
                        F2 =  fieldnames(S0.(F1{ifn_sub}));
                        for iifn = 1:length(F2(:));
                            if ~isfield(S(i).(F1{ifn_sub}),F2{iifn}) | isempty(S(i).(F1{ifn_sub})(1).(F2{iifn}))
                               S(i).(F1{ifn_sub})(1).(F2{iifn}) = S0.(F1{ifn_sub}).(F2{iifn});
                            end
                        end
                    end
                end
           end
           
           if ihp == 1; PP{it} = S;  end;
           if ihp == 2; HH{it} = S;  end;
       end
       
       
       %Field Content Check
       switch lower(TypePara{it})
%--------------------------------------------------------------------------------------------------
            case 'play'

                  P = PP{it}; n_P = length(P(:));
                  H = HH{it};

                  for i = 1:n_P
                      %NIO1 
                      nio.data = []; nio.sr = [];
                      P(i).nio1.type = lower(P(i).nio1.type);
                      if ~isempty(P(i).nio1.data); 
                         if isstr(P(i).nio1.data)
                            [j1, j2, FileExt] = fileparts(P(i).nio1.data);
                            switch lower(FileExt)
                                  case '.wav'
                                       try 
                                           [nio.data, nio.sr]= wavread(P(i).nio1.data); 
                                       catch; 
                                           P(i).nio1.status = -1; warning(sprintf('NIO1: Index %3.0d <<<< File: %s does not exist!',i, P(i).nio1.data)); 
                                           drawnow;
                                       end;
                                       
                                  otherwise
                                       try 
                                           nio = load(P(i).nio1.data); 
                                       catch; 
                                           H.nio1.status = -1; warning(sprintf('NIO1: Index %3.0d <<<< File: %s does not exist!',i, P(i).nio1.data));
                                           drawnow;
                                       end;
                            end
                            if ~isfield(nio,'sr') & isfield(nio,'t'); nio.sr = 1/(nio.t(2)-nio.t(1));end;
                            if ~isfield(nio,'data') | ~isfield(nio,'sr'); 
                                H.nio1.status = -1; warning(sprintf('NIO1: Index %3.0d <<<< Wrong File Formate: %s !',i, P(i).nio1.data));
                            else
                                P(i).nio1.sr = nio.sr;
                            end;
                         else
                            nio.data = P(i).nio1.data;
                         end
                      end
                      P(i).nio1.npt = size(nio.data,1);
                      if isempty(P(i).nio1.chan); 
                         P(i).nio1.chan = [0:size(nio.data,2)-1];
                      end
                      ind_chan = find(P(i).nio1.chan(:) >= 0);
                      P(i).nio1.nch = length(ind_chan(:));
                      if P(i).nio1.nch == 0 | P(i).nio1.npt == 0; 
                         warning(sprintf('NIO1: Index %3.0d <<<< No Data!',i)); 
                         drawnow;
                         P(i).nio1.status = 0; 
                      else;
                         P(i).nio1.status = 1;
                         if P(i).nio1.nch > size(nio.data,2); warning(sprintf('NIO1: Index %3.0d <<<< Missing Data Channel!',i));drawnow; H.nio1.status = -1; end;
                         if isempty(P(i).nio1.sr); warning(sprintf('NIO1: Index %3.0d <<<< No Sampling Rate!',i)); drawnow;  H.nio1.status = -1; end; 
                      end
                      clear('nio');
                      
                      %STG1
                      stg.data = []; 
                      if ~isempty(P(i).stg1.data); 
                         if isstr(P(i).stg1.data)
                            [j1, j2, FileExt] = fileparts(P(i).stg1.data);
                            switch FileExt
                                   case '.stg'
                                   otherwise
                                        try 
                                            stg = load(P(i).stg1.data); 
                                        catch; H.stg1.status = -1; warning(sprintf('STG1: Index %3.0d <<<< File: %s does not exist!',i, P(i).stg1.data)); drawnow; end;
                            end
                            if isfield(stg, 'sr'); P(i).stg1.sr = stg.sr; end;
                            if ~isfield(stg,'data') 
                                H.stg1.status = -1; warning(sprintf('STG1: Index %3.0d <<<< Wrong File Formate: %s !',i, P(i).stg1.data)); drawnow;
                            end;
                         else
                            stg.data  = P(i).stg1.data;
                         end
                      end
                      if isempty(P(i).stg1.chan); 
                         P(i).stg1.chan = [0:size(stg.data,2)-1];
                      end
                      for ich = 1:size(stg.data,2);
                          P(i).stg1.npt(ich) = length(stg.data{1,ich}(:));
                      end
                      ind_chan = find(P(i).stg1.chan >= 0); 
                      P(i).stg1.nch = length(ind_chan(:)); 
                      if P(i).stg1.nch == 0 | all(P(i).stg1.npt == 0); 
                         %warning(sprintf('STG1: Index %3.0d <<<< No Data!',i)); 
                         P(i).stg1.status = 0; 
                      else;
                         P(i).stg1.status = 1;
                         for ich = 1:size(stg.data,2)
                             if P(i).stg1.npt(ich) ~= length(stg.data{2,ich}(:)); 
                                 warning('STG1: Index %3.0d <<<< Time Stamp Error!'); 
                                 drawnow;
                                 P(i).stg1.status = -1; 
                              end;
                              if P(i).stg1.npt(ich) == 0 & P(i).stg1.chan(ich) >= 0; 
                                 warning('STG1: Index %3.0d <<<< Missing Data!'); 
                                 drawnow;
                                 P(i).stg1.status = -1; 
                              end;
                         end
                      end

                      %LPT1
                      if any(P(i).lpt1.code ~= 0);
                         P(i).lpt1.code(end+1:4) = 0; 
                         P(i).lpt1.status = 1; 
                      end
                      
                      %init P info text
                      if isempty(P(i).info)
                         info_file = '';
                         if isstr(P(i).nio1.data);  info_file = [info_file, P(i).nio1.data,' ;']; end;
                         if isstr(P(i).stg1.data);  info_file = [info_file, P(i).stg1.data,' ;']; end;
                         P(i).info = info_file;
                      end

                end
                
                nio1 = [P(:).nio1]; nid1 = [P(:).nid1]; lpt1 = [P(:).lpt1]; stg1 = [P(:).stg1]; com1 = [P(:).com1];
                if any([nio1(:).status] == 1);    H.nio1.status = 1; end;
                if any([nio1(:).status] == -1);   H.nio1.status = -1; end;
                if any([nid1(:).status] == 1);    H.nid1.status = 1; end;
                if any([nid1(:).status] == -1);   H.nid1.status = -1; end;
                if ~isempty(H.code_header);  H.lpt1.status = 1; end;
                if any([lpt1(:).status] == 1);    H.lpt1.status = 1; end;
                if any([lpt1(:).status] == -1);   H.lpt1.status = -1; end;
                if any([stg1(:).status] == 1);    H.stg1.status = 1; end;
                if any([stg1(:).status] == -1);   H.stg1.status = -1; end;
                if any([com1(:).status] == 1);    H.com1.status = 1; end;
                if any([com1(:).status] == -1);   H.com1.status = -1; end;

                %Init Header
                %---------------------------

                 H.type      = lower(H.type);
                 H.type_play = lower(H.type_play);
                     
                
                %Indices
                if isempty(H.SubInd); H.SubInd = [1:n_P]; end;
                [n, n_sub] = size(H.SubInd); size_sub = zeros(1,n_sub);
                for is = 1:n_sub; size_sub(1,is) = max(H.SubInd(:,is)); end;
                if size(H.SubInd,2) > 1
                   for i = 1:n; sub = num2cell(H.SubInd(i,:)); H.Ind(i) = sub2ind(size_sub,sub{:}); end;
                else
                   H.Ind    = H.SubInd; H.SubInd = [H.SubInd, ones(n,1)];
                end

                %Isi
                if isequal(H.isi_fix, inf);
                   H.Isi = ones(n,1)*inf;
                else
                   if length(H.isi_fix(:)) < n 
                      H.Isi = unifrnd(H.isi_fix(1), H.isi_fix+H.isi_var, n, 1);
                   elseif length(H.isi_fix(:)) == n
                      H.Isi = H.isi_fix;
                   end
                end

                %Playlist
                if length(H.list_play(:)) < n | ~iscell(H.list_play);
                   Count_stimrep = zeros(n_P,1);
                   H.list_play = cell([n,1]);
                   for i = 1:n; 
                       ind = H.Ind(i);
                       p = P(ind);
                       Count_stimrep(ind) = Count_stimrep(ind) + 1;
                       info_subind = ['( ', sprintf('%04u|', [H.SubInd(i,1:end-1)]), sprintf('%04u)', H.SubInd(i,end))];
                       H.list_play{i,1} = sprintf('(%04u|%04u|%04u); Sub-Index: %s; Info: %s',i, Count_stimrep(ind), ind-1, info_subind, p.info);
                   end
                end

                PP{it}   = P;
                HH{it}   = H;
%--------------------------------------------------------------------------------------------------
           case 'cftse'
                Cond = struct('type',[],'ln',[],'s1',[]);
                var_cond = {'type','ln','s1'};
                P = PP{it}; n_P = length(P(:));
                for ip = 1:n_P
                   P(ip).ncond = length(P(ip).s1(:));
                   P(ip).Cond = Cond;
                   P(ip).COND = {};
                   for ic = 1:P(ip).ncond
                       s1_str = P(ip).s1{ic};
                       P(ip).Cond(ic).type = P(ip).type;
                       P(ip).Cond(ic).s1 = s1_str;
                       
                       ind0 = findstr(s1_str,'LN(');
                       P(ip).Cond(ic).ln = strread(s1_str(ind0:end),['LN(%f'],1);
%                        ind0 = findstr(s1_str,'ID(');
%                        P(ip).Cond(ic).id = strread(s1_str(ind0:end),['ID(%f'],1);
%                       
%                        ind0 = findstr(s1_str,'US(');
%                        P(ip).Cond(ic).contingency = 0; 
%                        for ifi0 = 1:length(ind0(:))
%                          [cond_us(ifi0)]  = strread(s1_str(ind0(ifi0):end),'US(%s',1,'delimiter',',');
%                      
%                          switch cond_us{ifi0}
%                              case 'S'
%                                    P(ip).Cond(ic).contingency(ifi0) = -1;
%                              case 'N'
%                                    P(ip).Cond(ic).contingency(ifi0) = 1;
%                              case 'O'
%                                    P(ip).Cond(ic).contingency(ifi0) = 0;
%                          end
%                          
%                        end
%                        if ~all([P(ip).Cond(ic).contingency] == P(ip).Cond(ic).contingency(end)); P(ip).Cond(ic).contingency = 0;
%                        else; P(ip).Cond(ic).contingency = P(ip).Cond(ic).contingency(end); end;

                       P(ip).var_cond = var_cond;
                       P(ip).COND(ic,1:2) = {P(ip).Cond(ic).type,P(ip).Cond(ic).s1};
                       
                       P(ip).nrep(ic) = length(find(P(ip).s2 == P(ip).Cond(ic).ln ));
                      
                   end

               end               
               PP{it}   = P;
               
                             
%--------------------------------------------------------------------------------------------------
           case 'tone1'
                NV = {'f','att','dur'}; 
                Cond = struct('type',[],'f',[],'att',[],'dur',[]);
                n_var = length(NV(:));
                PV = cell([n_var,1]);
   
                P = PP{it}; n_P = length(P(:));
                for ip = 1:n_P
                   P(ip).chan_sig = P0.chan_sig;
                   P(ip).range_sig = P0.range_sig;
                   P(ip).type_trig = P0.type_trig;
                   P(ip).chan_trig = P0.chan_trig;
                   for iv = 1:n_var; 
                       PV{iv} = P(ip).(NV{iv});
                       N_var(iv) = length(PV{iv});
                       I_var{iv} = [1:N_var(iv)]';
                   end;
                   [I1, I2, I3] = ndgrid(I_var{:});
                   I = [I1(:), I2(:), I3(:)];
                   P(ip).ncond = size(I,1);
                   P(ip).COND = {};
                   P(ip).COND = repmat({P(ip).type},[P(ip).ncond,1]);
                   if P(ip).ncond > 0
                      for iv = 1:n_var; P(ip).COND(1:P(ip).ncond,iv+1) = num2cell(PV{iv}(I(:,iv))); end;
                   end
                   P(ip).var_cond = NV; 
                   P(ip).Cond = Cond;
                   for ic = 1:P(ip).ncond
                       P(ip).Cond(ic).type = P(ip).type;
                       for iv = 1:n_var; P(ip).Cond(ic).(NV{iv}) = P(ip).COND{ic,iv+1}; end;
                   end
               end               
               PP{it}   = P;
               
               
%--------------------------------------------------------------------------------------------------
           case 'puls1'

                NV = {'a_ph1','chan','dur_ph1'}; 
                Cond = struct('type',[],'a_ph1',[],'chan',[],'dur_ph1',[]);  
                n_var = length(NV(:));
                PV = cell([n_var,1]);
   
                P = PP{it}; n_P = length(P(:));
                for ip = 1:n_P
                   P(ip).type_trig = P0.type_trig;
                   P(ip).chan_trig = P0.chan_trig;
                   for iv = 1:n_var; 
                       PV{iv} = P(ip).(NV{iv});
                       N_var(iv) = length(PV{iv});
                       I_var{iv} = [1:N_var(iv)]';
                   end;
                   [I1, I2, I3] = ndgrid(I_var{:});
                   I = [I1(:), I2(:), I3(:)];
                   P(ip).ncond = size(I,1);
                   P(ip).COND = {};
                   P(ip).COND = repmat({P(ip).type},[P(ip).ncond,1]);
                   if P(ip).ncond > 0
                       for iv = 1:n_var; P(ip).COND(1:P(ip).ncond,iv+1) = num2cell(PV{iv}(I(:,iv))); end;
                   end
                   P(ip).var_cond = NV;
                   P(ip).Cond = Cond;
                   for ic = 1:P(ip).ncond
                       P(ip).Cond(ic).type = P(ip).type;
                       for iv = 1:n_var; P(ip).Cond(ic).(NV{iv}) = P(ip).COND{ic,iv+1}; end;
                   end
               end               
               PP{it}   = P;
%--------------------------------------------------------------------------------------------------
           case 'fmsweep'
                Cond = struct('type',[]);
                P = PP{it}; n_P = length(P(:));
                for ip = 1:n_P
                   P(ip).chan_sig  = P0.chan_sig;
                   P(ip).range_sig = P0.range_sig;
                   P(ip).type_trig = P0.type_trig;
                   P(ip).chan_trig = P0.chan_trig;
                   P(ip).ncond = 1;
                   P(ip).COND = {};
                   P(ip).COND(1:P(ip).ncond) = {P(ip).type};
                   P(ip).var_cond = {'type'};
                   P(ip).Cond = Cond; 
                   for ic = 1:P(ip).ncond
                       P(ip).Cond(ic).type = P(ip).type;
                   end
               end               
               PP{it}   = P;
%--------------------------------------------------------------------------------------------------
           otherwise
       end
       
    else
        keyboard
        HH{it} = H0; PP{it} = P0; FF0_P{it} = F0_P; FF0_H{it} = F0_H;
    end
    
end;

%output
if flag_def == 0;
   count_out = 0;
   for it = 1:count_type
       count_out = count_out + 1;
       varargout{count_out} = PP{it}; 
       count_out = count_out + 1;
       varargout{count_out} = FF0_P{it}; 
       if isequal(lower(TypePara{it}),'play')
          count_out = count_out + 1;
          varargout{count_out} = HH{it}; 
          count_out = count_out + 1;
          varargout{count_out} = FF0_H{it}; 
       end
   end
   
elseif flag_def == 1;
   count_out = 0;
   for it = 1:count_type
       count_out = count_out + 1;
       varargout{count_out} = PP{it}; 
       if isequal(lower(TypePara{it}),'play')
          count_out = count_out + 1;
          varargout{count_out} = HH{it}; 
       end
   end
end

