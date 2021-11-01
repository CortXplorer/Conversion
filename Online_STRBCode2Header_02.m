function [Header, Ind_header] = Online_STRBCode2Header_02(code, offset_code, flag_lh)

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
    n_field = (length(ind_field(:))-1)/2;
    
    % first create a cell array that contains the field names
    for i1 = 1:n_field
        code_Fi = code_H(ind_field((i1-1)*2+1)+1:ind_field((i1-1)*2+2)-1);
        ind_word = find(code_Fi == cw_word);
        n_word   = length(ind_word(:))-1;
        code_str = [];
        for i3 = 1 : n_word
            ci = num2str(code_Fi(ind_word(i3)+1:ind_word(i3+1)-1)-offset_code);
            is = find(isspace(ci)); ci(is) = [];
            code_str(i3) = str2num(ci(:)');
        end
            Fields{i1,1} = eval(char(code_str));
    end

    % now get the contents of the fields
    for i2 = 1:n_field
        code_Fi = code_H(ind_field((i2-1)*2+2)+1:ind_field((i2-1)*2+3)-1);
        ind_word = find(code_Fi == cw_word);
        n_word   = length(ind_word(:))-1;
        code_str = [];
        for i3 = 1 : n_word
            ci = num2str(code_Fi(ind_word(i3)+1:ind_word(i3+1)-1)-offset_code);
            is = find(isspace(ci)); ci(is) = [];
            code_str(i3) = str2num(ci(:)');
        end
        if ~isempty(code_str)
            try
            Header(i).(Fields{i2,1}) = eval(char(code_str));
            catch
                try % here is a nasty hack b/c some cellmix arrays have been treated as cellstr
                    test = char(code_str);
                    i_notstr = find(code_str < 32 | code_str > 127);
                    for i1 = 1:size(i_notstr,2)
                        test = [test(1:i_notstr(end-i1+1)-1),num2str(code_str(i_notstr(end-i1+1))),test(i_notstr(end-i1+1)+1:end)];
                    end
                    Header(i).(Fields{i2,1}) = eval(test);
                catch
                    keyboard;
                end
            end   
        else
            Header(i).(Fields{i2,1}) = [];
        end
    end
end

if flag_lh == 1 & n_header > 0; Header = Header(end); end;