function [Header, Ind_header] = Play_Code2Header_01(code, offset_code, flag_lh)

if ~exist('offset_code') | isempty(offset_code); offset_code  = 255-16; end;
if ~exist('flag_lh') | isempty(flag_lh); flag_lh  = 0; end;
cw_header  = offset_code+14;
cw_field   = offset_code+13;
cw_word    = offset_code+12;

%init
Header = [];
Fields = {'type','string'};
Ind_header = [];

%find header 
ind_header = find(code(:) == cw_header);
n_header = floor(length(ind_header(:))/2);

for i = 1:2:n_header*2;
    code_H = code(ind_header(i)+1:ind_header(i+1)-1);
    Ind_header = [Ind_header; [ind_header(i):ind_header(i+1)]'];
    
    %get type
    ind_field = find(code_H == cw_field);
    n_field = length(ind_field(:));
    for ii = 1:n_field-1
        code_Fi = code_H(ind_field(ii)+1:ind_field(ii+1)-1);
        ind_word = find(code_Fi == cw_word);
        n_word   = length(ind_word(:))-1;
        code_str = [];
        for iii = 1 : n_word
            ci = num2str(code_Fi(ind_word(iii)+1:ind_word(iii+1)-1)-offset_code);
            is = find(isspace(ci)); ci(is) = [];
            code_str(iii) = str2num(ci(:)');
        end
        if ~isempty(code_str)
            Header(i).(Fields{ii,1}) = eval(char(code_str));
        else
            Header(i).(Fields{ii,1}) = [];
        end
        if ii == 1
           [j, Fields] = Play_ParaCheck_02(Header.type); 
        end
    end
    [Header(i)] = Play_ParaCheck_02(Header.type, Header); 
end

if flag_lh == 1 & n_header > 0; Header = Header(end); end;




