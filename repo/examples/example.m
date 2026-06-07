% Minimal usage example for the DNAm Network Clock (MATLAB).
%
% Run from the repository root:
%
%     >> example
%
% Or copy these commands into your own script.

addpath('matlab');

% Load reference input shipped with the repo
T = readtable(fullfile('examples', 'example_input.csv'), ...
              'ReadRowNames', true, 'VariableNamingRule', 'preserve');
cpg_ids = T.Properties.VariableNames;
betas   = table2array(T);
fprintf('Loaded reference input: %d samples × %d CpGs\n', ...
        size(betas, 1), size(betas, 2));

% Predict
ages = predict_age(betas, cpg_ids);
fprintf('\nPredictions:\n');
fprintf('  range  : %.2f – %.2f years\n', min(ages), max(ages));
fprintf('  median : %.2f years\n', median(ages));

% Verify against stored expected output (Python-generated)
expected = readtable(fullfile('examples', 'expected_output.csv'));
max_diff = max(abs(ages - expected.predicted_age));
if max_diff < 1e-6
    fprintf('\n  max |MATLAB − expected| = %.2e   ✓ PASS\n', max_diff);
else
    fprintf('\n  max |MATLAB − expected| = %.2e   ✗ FAIL\n', max_diff);
end
