function data = read_csvfile(filename)
% READ_CSVFILE reads a CSV file from the PPG/ECG CSV files.

% determine CSV directory relative to this script.
cwd = fileparts(mfilename('fullpath'));
dir = fullfile(cwd, 'bidmc-ppg-and-respiration-dataset-1.0.0', 'bidmc_csv');
if (exist(dir, 'dir') ~= 7)
    error('CSV directory missing!');
end

% simplify the name parsing.
filename = strrep(filename, '_Signals.csv', '');
filename = [filename, '_Signals.csv'];
sigfile = fullfile(dir, filename);
if (exist(sigfile, 'file') ~= 2)
    error('Specified FILENAME "%s" is not a valid CSV.');
end
data = readtable(sigfile);

end