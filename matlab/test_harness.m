function test_harness(results_filename)
% TEST_HARNESS will kick off the filtering test harness using Monte-Carlo
% simulations.
%
%   TEST_HARNESS(RESULTS_FILENAME)
%       RESULTS_FILENAME (char) is the output .MAT file to save the results
%       to. If an existing filename is found, we'll increment the name
%       counter so that no old results are saved over.
%
% See Also: ACN, MAWGN, EMD_RILLING, WNOISE
%
% Contact: Jordan R. Smith

if (nargin < 1)
    % "results_10Apr2022_1026.mat" for example.
    results_filename = sprintf('results_%s.mat', ...
        strrep(datestr(now, 'ddmmmyy HHMMSS'), ' ', '_'));
end
workdir = fileparts(mfilename('fullpath'));
outfile = fullfile(workdir, results_filename);


testenv = struct;                       % primary test container.
testenv.params.snr = -12:2:8;           % SNR contamination axis
testenv.params.z = 0.95;                % Confidence Interval z-score
testenv.params.use_parallel = true;     % use parallel-processing to speed up?
testenv.params.mintrials = 20;          % minimum # of trials before checking test conditions
testenv.params.noisetype = 'pink';      % noise type to contaminate with
testenv.params.N = 14;                  % 2^N signal length
testenv.params.tolerance = 0.5;         % error tolerance in Monte-Carlo loop
testenv.params.date = datestr(now, ...  % encoding time of the test.
    'dd-mmm-yyyy HH:MM:SS');

% Loop each signal and test.
testenv.signals = struct('input', {}, 'label', {}, 'results', {});

testenv.signals(1).input = wnoise(3, testenv.N);
testenv.signals(1).label = 'Heavy sine';
testenv.signals(1).results = struct('snr', [], 'mse', []);

testenv.signals(2).input = wnoise(4, testenv.N);
testenv.signals(2).label = 'Doppler';
testenv.signals(2).results = struct('snr', [], 'mse', []);

% save the results.
save(outfile, 'testenv');


end % test_harness

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_routine(env, signal_index)
% TEST_ROUTINE is the test loop container.
%
%   TEST_ROUTINE(ENV, SIGNAL_INDEX)
%       ENV (struct) is the test container.
%       SIGNAL_INDEX (scalar) indicates which signal in the struct array to
%           gather results for.
%
% See Also: ACN

    % splat variables from "env" struct. If parallel-processing is used to
    % speed up execution, we don't want the compiler to memcpy the entire
    % container.
    ntype = env.params.noisetype;
    xin = env.params.signal(signal_index).input;
    snr = env.params.snr;

    if (env.params.use_parallel)
        % parallel-processing, speeding up with CPU core workers.
        parfor (pIdx = 1:length(env.params.snr))
            iMc = 1;
            while(1)
                switch ntype
                    case {'pink', 'violet', 'red', 'blue'}
                        xnoise = ACN(xin, snr(pIdx), ntype);
                    case 'white'
                        xnoise = mAWGN(xin, snr(pIdx));
                    otherwise
                        error('Unknown noise type "%s".', env.params.noisetype);
                end % switch
                
                % TODO: perform the filtering and MC check.

            end % while
        end % par-for
    else
        % no parallel-processing, conventional for-loop testing.
        for idx = 1:length(env.params.snr)
            iMc = 1;
            while(1)
                switch env.params.noisetype
                    case {'pink', 'violet', 'red', 'blue'}
                        xnoise = ACN(xin, snr(idx), ntype);
                    case 'white'
                        xnoise = mAWGN(xin, snr(idx));
                    otherwise
                        error('Unknown noise type "%s".', env.params.noisetype);
                end % switch

                % TODO: perform the filtering and MC check.

            end % while

        end % for
    end % if

end % test_routine

                




                    

