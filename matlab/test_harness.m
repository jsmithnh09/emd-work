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
partialdir = fullfile(workdir, 'partial');
if (exist(partialdir, 'dir') ~= 7)
    mkdir(partialdir);
end


testenv = struct;                       % primary test container.
testenv.params.snr = -12:2:8;           % SNR contamination axis
testenv.params.z = 0.95;                % Confidence Interval z-score
testenv.params.mintrials = 20;          % minimum # of trials before checking test conditions
testenv.params.noisetype = 'pink';      % noise type to contaminate with
testenv.params.N = 14;                  % 2^N signal length
testenv.params.tolerance = 0.5;         % error tolerance in Monte-Carlo loop
testenv.params.date = datestr(now, ...  % encoding time of the test.
    'dd-mmm-yyyy HH:MM:SS');


% Loop each signal and test.
testenv.signals = struct('input', {}, 'label', {}, 'filters', {});

% pre-allocate results based on length of SNR axis.
res_axis = repmat(-Inf, 1, length(testenv.params.snr));
res_struct = struct('snr', res_axis, 'mse', res_axis);

filt_struct = struct('func', {@emdB, @emdK, @emdH, @emdHD}, ...
    'label', {'EMD-Beta', 'EMD-Kurtosis', 'EMD-H', 'EMD-Hausdorff'}, ...
    'result', {res_struct, res_struct, res_struct, res_struct});


% declare the signals that ought to be tested.
testenv.signals(1).input = wnoise(3, testenv.params.N);
testenv.signals(1).label = 'Heavy sine';
testenv.signals(1).filters = filt_struct;

testenv.signals(2).input = wnoise(4, testenv.params.N);
testenv.signals(2).label = 'Doppler';
testenv.signals(2).filters = filt_struct;

% pull in an ECG signal
[t, ecg, ~] = read_csvfile('bidmc_01');
fs = 1/(t(2)-t(1));
start = 3 * fs; % determine start point (arbitrary... 3 seconds?)
stop = start + (2^testenv.params.N); % keep it the same length as the other test signals.
ecg = ecg(start:stop-1);

testenv.signals(3).input = ecg;
testenv.signals(3).label = 'ECG';
testenv.signals(3).filters = filt_struct;

fprintf(1, '%s\n', repmat('*', 40, 1));
fprintf(1, ' Test started: %s \n', datestr(now, 'dd-mmm-yyyy HH:MM PM'));
for testIdx = 1:length(testenv.signals)

    % indicate which signal.
    fprintf(1, 'Testing %s\n', testenv.signals(testIdx).label);
    testenv = test_routine(testenv, testIdx);

    % because each of these signals can take a considerable amount of time,
    % lets save out a partial file in case we need to stop at any point.
    fname = sprintf('partial_signal%d_%s.mat', testIdx, ...
        datestr(now, 'mm-dd-HHMMSS'));
    save(fullfile(partialdir, fname), 'testenv', '-mat');


end


% save the results.
save(outfile, 'testenv');


end % test_harness

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function env = test_routine(env, signal_idx)
% TEST_ROUTINE is the test loop container.
%
%   ENV = TEST_ROUTINE(ENV, SIGNAL_INDEX)
%       ENV (struct) is the test container.
%       SIGNAL_INDEX (scalar) indicates which signal in the struct array to
%           gather results for.
%
% The test routine utilizes Monte Carlo testing across the SNR axis stored
% in the environment struct.
%
% See Also: ACN

% splat variables from "env" struct. If parallel-processing is used to
% speed up execution, we don't want the compiler to memcpy the entire
% container.
ntype = env.params.noisetype;
xin = env.signals(signal_idx).input;
snr = env.params.snr;
tol = env.params.tolerance;
mintrials = env.params.mintrials;

% no parallel-processing, conventional for-loop testing.
for snr_idx = 1:length(env.params.snr)
    try
        fprintf(1, 'Testing SNR: %.3g...', env.params.snr(snr_idx))
        for filt_idx = 1:length(env.signals(signal_idx).filters)
            filter_func = env.signals(signal_idx).filters(filt_idx).func;
            filter_name = env.signals(signal_idx).filters(filt_idx).label;
            mcIdx = 1;
            [mcSNR, mcMSE, sigma, ci] = deal([]);

            while(1)

                % generate/apply noise with SNRdB and noise-type.
                xn = acn(xin, snr(snr_idx), ntype);

                % extract the IMF and perform de-noising.
                ximf = emd_rilling(xn);
                try
                    yout = filter_func(ximf);
                catch ME
                    error('test_harness:routine:: "%s" filter error SNR=%g Monte Carlo #(%d):\n%s', ...
                        filter_name, snr(snr_idx), mcIdx, ME.message);
                end

                % tracking performance/CI based on SNR.
                [mcSNR(mcIdx), mcMSE(mcIdx)] = perf(xin, yout);
                sigma(mcIdx) = std(mcSNR);

                % check if confidence limit exceeds interval.
                ci(mcIdx) = 1 - 2*qfunc(sqrt(mcIdx)/sigma(mcIdx)*tol);
                if (ci(mcIdx) >= 0.95) && (mcIdx > mintrials)
                    break;
                else
                    mcIdx = mcIdx + 1;
                    continue;
                end % if
            end % while

            % store the results.
            env.signals(signal_idx).filters(filt_idx).result.snr(snr_idx) = mean(mcSNR);
            env.signals(signal_idx).filters(filt_idx).result.mse(snr_idx) = mean(mcMSE);

        end % filter-for
    catch ME
        fprintf(1, 'ERROR\n');
        rethrow(ME);
    end
    fprintf(1, 'done\n');
end % SNR-for

end % test_routine








