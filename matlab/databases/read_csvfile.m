function [t, ppg, ecg] = read_csvfile(filename)
% READ_CSVFILE reads a CSV file from the PPG/ECG CSV files.
%
%   [T, PPG, ECG] = READ_CSV(FILENAME)
%       FILENAME (char) indicates the filename pattern. The files  are
%       annoted "bidmc_%d_Signals.csv", where %d denotes a number.
%
% Example:
%   [t, ppg, ecg] = read_csvfile('bidmc_01');
%   Fs = 1 / (t(2)-t(1)); % 125 Hz.
%
% See Also: READTABLE

% determine CSV directory relative to this script.
name = 'bidmc-ppg-and-respiration-dataset-1.0.0';
cwd = fileparts(mfilename('fullpath'));
dir = fullfile(cwd, name);
if (exist(dir, 'dir') ~= 7) && ...
        (exist(fullfile(cwd, [name, '.zip']), 'file') == 2)
    % indicates the ZIP is available, but not unzipped yet...
    unzip(fullfile(cwd, [name, '.zip']));
end

if (exist(dir, 'dir') ~= 7)
    error('CSV directory missing or the ZIP file hasn''t been unzipped.');
end

% simplify the name parsing.
filename = strrep(filename, '_Signals.csv', '');
filename = [filename, '_Signals.csv'];
sigfile = fullfile(dir, filename);
if (exist(sigfile, 'file') ~= 2)
    error('Specified FILENAME "%s" is not a valid CSV.');
end

% disable the warning state for renaming variables.
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
dtable = readtable(sigfile);
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');
t = dtable.Time_s_;
ppg = dtable.PLETH;
ecg = dtable.II;
clear('dtable');
return;

end