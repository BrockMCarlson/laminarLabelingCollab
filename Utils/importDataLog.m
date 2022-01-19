function TuneList = importTuneList(flag_dionly)
%   case {1,'all di tasks','dionly','ditasks'}
%   case {2,'brfs','Brock'} % Brock
%   case {3,'cosinteroc','mcosinteroc','contrast','Blake'} % Brock % Blake
%   case {4,'bmcBRFS'};
if nargin < 1
    flag_dionly = 0; 
end

%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /ephys-analysis/stim/DataLog.xlsx
%    Worksheet: TuneList

%% Import the data
if exist('DataLog.xlsx','file')
    [~, ~, raw] = xlsread('DataLog.xlsx','TuneList');
elseif exist('./stim/DataLog.xlsx','file')
     [~, ~, raw] = xlsread('/stim/DataLog.xlsx','TuneList');
elseif exist('./ephys-analysis/stim/DataLog.xlsx','file')
     [~, ~, raw] = xlsread('DataLog.xlsx','TuneList');
else
    global HOMEDIR
    if ~isempty(HOMEDIR) && exist( [HOMEDIR '/ephys-analysis/stim/DataLog.xlsx'],'file')
        [~, ~, raw] = xlsread([HOMEDIR '/ephys-analysis/stim/DataLog.xlsx'],'TuneList');
        clear HOME
    else
        clear global HOME
        error('cannot find DataLog.xls');
    end
end
columlim = 29;
headersrow = 1;
headers =  raw(headersrow,1:columlim);
raw = raw(headersrow+1:end,1:columlim);
rowlim = find(cellfun(@isnan,raw(:,1)),1) - 1; 
raw = raw(1:rowlim,:); 

%% Allocate imported array to headers variable names

for f = 1:length(headers)
    dat = raw(:,f);
    
    switch headers{f}
        case 'Date'
            dat = cellfun(@num2str,dat,'UniformOutput',0);
            TuneList.('Datestr') = dat;
            TuneList.('Datenum') = cell2mat(cellfun(@(x) datenum(x,'yymmdd'),dat,'UniformOutput',0));
        case {'SinkBtm','Rig','Drobo','BadBtmCh','GridX','GridY','Rank'}
            I = cellfun(@(x) ~isnumeric(x) && ~islogical(x),dat);
            dat(I) = {NaN};
            dat = cell2mat(dat);
            TuneList.(headers{f}) = dat;
        case {'dotmapping','rfori','rfsf','rfsize','drfori','cosinteroc','mcosinteroc','dmcosinteroc','brfs','dbrfs','rsvp','evp','darkrest','rfsfdrft','bminteroc','bmcBRFS'}
            TuneList.(headers{f}) = cellfun(@eval,dat,'UniformOutput',0);
        otherwise
            TuneList.(headers{f}) = dat;
    end
end

% remove bad penetrations based on rank
exclude = TuneList.Rank == 0;
fields = fieldnames(TuneList);
for f = 1:length(fields)
    TuneList.(fields{f})(exclude) = [];
end

% if instructed, remove sessions w/o di tasks
if isstring(flag_dionly) || (isnumeric(flag_dionly) && flag_dionly > 0)
    clear I
    switch flag_dionly
        case {1,'all di tasks','dionly','ditasks'}
            ditasks = {'cosinteroc','mcosinteroc','dmcosinteroc','brfs','dbrfs', 'rsvp'};
        case {2,'brfs','Brock'} % Brock
            ditasks = {'brfs'};
        case {3,'blake'} % Brock % Blake
            ditasks = {'cosinteroc','mcosinteroc','brfs'};
        case {4,'bmcBRFS'} % Neural data recorded by BMC and BM, 2021
            ditasks = {'bmcBRFS'}; % This is just making sure the session has ____ before analyzing all paradigms.
    end
    
    for d = 1:length(ditasks)
        I(:,d) = cellfun(@isempty,TuneList.(ditasks{d}));
    end
    I = all(I,2);
    
    fields = fieldnames(TuneList);
    for f = 1:length(fields)
        TuneList.(fields{f})(I) = [];
    end
end



