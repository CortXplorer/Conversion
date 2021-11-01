function DS  = Convert_Header2SWEEP_01(Header,t_pre,dur,stim_list)

% ------------ check for experiment type --------------------------------
switch Header.exp_type 
    case 'Single Tones'
        % use structure of raw wave files coding
        DS.type        = Header.exp_type;
        DS.t_Pre       = t_pre;
	    DS.t_Sig       = Header.stimlist(:,end-1)';
	    DS.t_Pau       = dur-t_pre;
	    DS.ISI         = Header.ISI;
        DS.Parm0       = Header.stimlist(:,end-1)';	% Attenuation
        DS.Parm1       = 1;	% min. stimulus index = 1
        DS.Parm2       = Header.ncond;	% max. stimulus index = DS.n_Stm
        DS.Parmd       = 1;	% index step = 1
        DS.n_Rep	   = Header.n_rep;	% number of repetitions
        DS.n_Stm	   = Header.ncond;
       
    otherwise % right now it seems that this universally usable so just copy the code
        % use structure of raw wave files coding
        DS.type        = Header.exp_type;
        DS.t_Pre       = t_pre;
	    DS.t_Sig       = Header.stimlist(:,end-1);
	    DS.t_Pau       = dur-t_pre;
	    DS.ISI         = Header.ISI;
        DS.Parm0       = Header.stimlist(:,end-1);	% Attenuation
        DS.Parm1       = 1;	% min. stimulus index = 1
        DS.Parm2       = Header.ncond;	% max. stimulus index = DS.n_Stm
        DS.Parmd       = 1;	% index step = 1
        DS.n_Rep	   = Header.n_rep;	% number of repetitions
        DS.n_Stm	   = Header.ncond;
end




