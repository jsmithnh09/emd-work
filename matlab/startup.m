function startup
%STARTUP file for initializing MATLAB on startup.
%
%   STARTUP will be called on starting MATLAB to add the appropriate files
%   to the $PATH for common use.
%
% Contact: Jordan R. Smith

[~] = fprintf('*** emd-work Configuration ***\n');
basedir = fileparts(mfilename('fullpath'));
tbxlist = {'common', 'emd', 'ggd'};

for iTbx = 1:length(tbxlist)
    tbxdir = fullfile(basedir, tbxlist{iTbx});
    if (exist(tbxdir, 'dir') ~= 7)
        warning('startup:MissingTbx', 'Missing toolbox "%s".', tbxdir);
    else
        addpath(tbxdir);
    end
end

end % startup

