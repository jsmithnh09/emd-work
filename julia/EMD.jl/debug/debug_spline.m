% pull the bad iteration in.

raw = fileread('../../src/debug_badinput.txt');
raw = strsplit(raw, newline);
raw = raw(1:end-1);
bad = cellfun(@str2double, raw);
t = 1:length(bad);
tol = [0.05, 0.5, 0.05];

% check the min/max
[indmin, indmax, indzer] = extr(bad,1:length(bad));

% get the spline inputs.
[tmin, tmax, mmin, mmax] = boundary_conditions(indmin, indmax, t, bad, bad, 2);

% perform the spline.
envmin = spline(tmax, mmax, t);