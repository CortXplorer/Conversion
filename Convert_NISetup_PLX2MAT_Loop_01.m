function Convert_NISetup_PLX2MAT_Loop_01

% Some useful variables
path_in_dat = 'D:\MyCode\Conversion';
path_default_plx = 'D:\MyCode\Conversion';

% Identify path
path_current = cd;
cd(path_default_plx)
[name_file, name_path] = uigetfile('*.plx', ...
    'Please select an arbitrary file to identify folder.');
cd(path_current);

% Reconstruct file names
files1 = dir([name_path, '*.plx']);
filenames = char({files1.name}');
n_files = size(filenames,1);

% loop only through files not yet analysed
files=dir([name_path, '*.plx']);
files3=dir([name_path, '*.mat']);
names = sort({files.name}');
names2=sort({files3.name}');
names2=strrep(names2, '.mat', '.plx');

x = [];
for n = 1:size(names2,1)
    x = [x;strmatch(names2(n,:), names, 'exact')];
end
loop = setxor(x,[1:size(names,1)]);
names = names(loop);
n_files = length(names);

answers=inputdlg({'t_pre [s]','dur [s]'},'Choose default conversion parameters',1,{'0.2','0.7'});
t_pre  = str2num(answers{1});
dur    = str2num(answers{2});

% Loop through the chosen dat and their corresponding plx-files
    for i = 1 : n_files
        name_in_plx = char(names(i,:));
        disp(        '*************************************************');
        fprintf('Processing file   : %d of %d\n', i, n_files);
        fprintf('Analyzing         : %s\n', name_in_plx);
        fprintf('in folder         : %s\n', name_path);
        disp(        '*************************************************');
        if (isempty(name_in_plx) ~= 1)           
             Convert_NISetup_PLX2MAT_04([name_path, name_in_plx],1,1,t_pre,dur);
        else 
            disp(' ');
            disp('No corresponding *.plx-file chosen!');
            disp(' ');
        end    
        
    end
