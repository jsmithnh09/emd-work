function defaults = optparse(defaults, varargin)
% OPTPARSE performs optional keyword argument parsing.
%
%  OPTPARSE(DEFAULTS, 'NAME1', VALUE1, ...) accepts a default struct with
%  all possible Name/Value pairs, and the VARARGIN component is then
%  iterated to overwrite any default parameters.
%
%  This would allow a top-level function to accept mutliple key-words and
%  use the fieldnames of the DEFAULTS struct to determine what is/isn't
%  acceptable.
%
% See Also: FIELDNAMES, ISMEMBER

if (~isstruct(defaults))
    error('optparse:InvalidValue', ...
        'Input DEFAULTS must be a default parameter struct.');
end
fnames = fieldnames(defaults);
nargs = numel(varargin);
if (mod(nargs, 2))
    error('optparse:InvalidSyntax', ...
        'Input NAME/VALUE must be in pairs. Odd number of arguments.');
end
keys = varargin(1:2:end-1);
vals = varargin(2:2:end);
if (~all(cellfun(@ischar, keys)))
    error('optparse:InvalidValue', ...
        'All keywords must be type CHAR.');
end
present = ismember(keys, fnames);
if any(~present)
    error('optparse:InvalidValue', ...
        'Unknown keyword(s) "%s".', strjoin(keys, ', '));
end
for iF = 1:numel(keys)
    defaults.(keys{iF}) = vals{iF};
end

end % optparse
