function ages = predict_age(betas, cpg_ids, weights_dir)
%PREDICT_AGE  Apply the K=60 DNAm Network Clock to a methylation beta matrix.
%
%   This is a single, self contained file. There is nothing to install. Put it
%   on your MATLAB path, or just keep it next to the "weights" folder, then
%   call it from your own script.
%
%   ages = PREDICT_AGE(betas, cpg_ids)
%   ages = PREDICT_AGE(betas, cpg_ids, weights_dir)
%
%   Inputs
%     betas       Numeric matrix of size [n_samples x n_cpgs]. Each ROW is one
%                 sample, each COLUMN is one CpG. NaN values are tolerated.
%     cpg_ids     Cell array (or string array) of length n_cpgs giving the CpG
%                 identifier (Illumina cg string) for each column of betas.
%                 Missing CpGs are tolerated.
%     weights_dir (Optional) Path to the "weights" folder. If omitted, the
%                 function looks for ../weights relative to its own location
%                 (the repository layout), then ./weights, then pwd/weights.
%
%   Output
%     ages        [n_samples x 1] column vector of predicted ages (years), in
%                 the same row order as betas.
%
%   Example
%     T = readtable('examples/example_input.csv', 'ReadRowNames', true, ...
%                   'VariableNamingRule', 'preserve');
%     cpgs  = T.Properties.VariableNames;
%     betas = table2array(T);
%     ages  = predict_age(betas, cpgs);
%
%   The clock is a K=60 principal component ridge regression on per module
%   medians of beta values. See README.md for the methodology, the trained
%   model metrics, and the paper reference.
%
%   Author: Anton Carcedo. Licence: MIT.

% Argument handling.
if nargin < 2
    error('predict_age:NotEnoughInputs', ...
        'Provide both a beta matrix and a cell array of CpG IDs.');
end
if nargin < 3 || isempty(weights_dir)
    here = fileparts(mfilename('fullpath'));
    candidates = {fullfile(here, '..', 'weights'), ...
                  fullfile(here, 'weights'), ...
                  fullfile(pwd, 'weights')};
    weights_dir = '';
    for k = 1:numel(candidates)
        if isfolder(candidates{k})
            weights_dir = candidates{k};
            break;
        end
    end
    if isempty(weights_dir)
        error('predict_age:WeightsNotFound', ...
            ['Could not locate the "weights" directory automatically. Pass ' ...
             'the path as the third argument.']);
    end
end

if size(betas, 2) ~= numel(cpg_ids)
    error('predict_age:DimMismatch', ...
        'numel(cpg_ids) (%d) must equal number of columns in betas (%d).', ...
        numel(cpg_ids), size(betas, 2));
end

% Load weights, build the module matrix, apply the linear model.
W = load_clock_weights(weights_dir);
X = build_module_matrix(betas, cpg_ids, W);     % [n_samples x n_modules]

X = bsxfun(@minus,   X, W.scaler_mean');         % standardise
X = bsxfun(@rdivide, X, W.scaler_scale');
X = bsxfun(@minus,   X, W.pca_mean');            % centre for PCA

PC_scores = X * W.pca_components';               % [n_samples x 60]
ages = PC_scores * W.pc_coefs + W.intercept;     % [n_samples x 1]
end


% =====================================================================
% Local function: load and cache the K=60 clock weights from the CSVs.
% =====================================================================
function W = load_clock_weights(weights_dir)
%   Returns a struct W in canonical module column order (length 4595 unless
%   otherwise stated):
%     .module_order        [4595x1 double]  module IDs
%     .module_to_cpgs      {4595x1 cell}    CpG IDs per module
%     .module_to_strategy  {4595x1 cell}    'median' or 'singleton'
%     .scaler_mean         [4595x1 double]
%     .scaler_scale        [4595x1 double]
%     .pca_mean            [4595x1 double]
%     .pca_components      [60x4595 double]  PC loadings (rows = PCs)
%     .pc_coefs            [60x1 double]     ridge coefficients
%     .intercept           scalar            ridge intercept
%     .K                   scalar = 60
%     .n_modules           scalar = 4595
%
%   Results are cached in a persistent variable so repeated calls do not
%   re-read the CSVs.

    persistent cache cache_key
    key = weights_dir;
    if ~isempty(cache) && strcmp(key, cache_key)
        W = cache; return;
    end

    % module_definitions.csv
    defs_path = fullfile(weights_dir, 'module_definitions.csv');
    if ~isfile(defs_path)
        error('load_clock_weights:NoFile', 'Missing %s', defs_path);
    end
    defs = readtable(defs_path, 'TextType', 'string');

    % module_params.csv
    params = readtable(fullfile(weights_dir, 'module_params.csv'));
    module_order = params.Module;
    n_modules    = numel(module_order);

    % Per module CpG lists, in module_order.
    module_to_cpgs     = cell(n_modules, 1);
    module_to_strategy = cell(n_modules, 1);
    [tf, loc] = ismember(defs.Module, module_order);
    if ~all(tf)
        error('load_clock_weights:UnknownModule', ...
            'module_definitions.csv contains modules not in module_params.csv.');
    end
    for j = 1:n_modules
        sel = (loc == j);
        module_to_cpgs{j}     = cellstr(defs.CpG(sel));
        s = unique(defs.Strategy(sel));
        module_to_strategy{j} = char(s(1));
    end

    % pca_components_k60.csv, drop the 'PC' label column -> [60 x 4595].
    comps_tbl = readtable(fullfile(weights_dir, 'pca_components_k60.csv'), ...
                          'VariableNamingRule', 'preserve');
    pca_components = table2array(comps_tbl(:, 2:end));
    if size(pca_components, 1) ~= 60 || size(pca_components, 2) ~= n_modules
        error('load_clock_weights:BadComponents', ...
            'pca_components_k60.csv has shape %dx%d, expected 60x%d.', ...
            size(pca_components, 1), size(pca_components, 2), n_modules);
    end

    % ridge_k60.csv
    ridge = readtable(fullfile(weights_dir, 'ridge_k60.csv'), 'TextType', 'string');
    pc_mask   = startsWith(ridge.name, "PC");
    pc_coefs  = ridge.value(pc_mask);
    intercept = ridge.value(ridge.name == "intercept");
    if numel(intercept) ~= 1
        error('load_clock_weights:BadRidge', ...
            'ridge_k60.csv must contain exactly one row with name="intercept".');
    end

    % Assemble and cache.
    W = struct(...
        'module_order',       module_order, ...
        'module_to_cpgs',     {module_to_cpgs}, ...
        'module_to_strategy', {module_to_strategy}, ...
        'scaler_mean',        params.scaler_mean, ...
        'scaler_scale',       params.scaler_scale, ...
        'pca_mean',           params.pca_mean, ...
        'pca_components',     pca_components, ...
        'pc_coefs',           pc_coefs, ...
        'intercept',          intercept, ...
        'K',                  60, ...
        'n_modules',          n_modules);
    cache = W;
    cache_key = key;
end


% =====================================================================
% Local function: aggregate a samples x CpGs matrix into samples x modules.
% =====================================================================
function X = build_module_matrix(betas, cpg_ids, W)
%   Each output column is the per module value: median of the available CpGs
%   for 'median' modules, the single CpG's value for 'singleton' modules, or
%   the training mean if a module has no input CpGs at all. Emits a warning if
%   more than 5% of clock CpGs are absent from cpg_ids.

    n_samples = size(betas, 1);
    n_modules = W.n_modules;
    X = zeros(n_samples, n_modules);

    % CpG ID -> column index in betas (case sensitive).
    if isstring(cpg_ids), cpg_ids = cellstr(cpg_ids); end
    cpg_to_col = containers.Map(cpg_ids, num2cell(1:numel(cpg_ids)));

    n_missing_modules    = 0;
    n_singleton_fallback = 0;
    n_cpgs_used          = 0;
    n_cpgs_needed        = 0;

    for j = 1:n_modules
        cpgs     = W.module_to_cpgs{j};
        strategy = W.module_to_strategy{j};
        n_cpgs_needed = n_cpgs_needed + numel(cpgs);

        % Resolve which CpGs are present in the input.
        present_cols = zeros(1, numel(cpgs));
        n_present = 0;
        for k = 1:numel(cpgs)
            if isKey(cpg_to_col, cpgs{k})
                n_present = n_present + 1;
                present_cols(n_present) = cpg_to_col(cpgs{k});
            end
        end
        present_cols = present_cols(1:n_present);
        n_cpgs_used = n_cpgs_used + n_present;

        if strcmp(strategy, 'singleton')
            if n_present >= 1
                X(:, j) = betas(:, present_cols(1));
            else
                X(:, j) = W.scaler_mean(j);
                n_singleton_fallback = n_singleton_fallback + 1;
            end
        else  % 'median'
            if n_present == 0
                X(:, j) = W.scaler_mean(j);
                n_missing_modules = n_missing_modules + 1;
            else
                if n_present == 1
                    X(:, j) = betas(:, present_cols(1));
                else
                    X(:, j) = median(betas(:, present_cols), 2, 'omitnan');
                end
                nan_rows = isnan(X(:, j));
                if any(nan_rows)
                    X(nan_rows, j) = W.scaler_mean(j);
                end
            end
        end
    end

    coverage = n_cpgs_used / n_cpgs_needed;
    if (1 - coverage) > 0.05
        warning('predict_age:MissingCpGs', ...
            ['DNAm Network Clock: %.1f%% of the clock''s CpGs are missing ' ...
             'from the input (%d/%d present). %d modules were filled with ' ...
             'the training mean. Predictions may be biased toward the ' ...
             'training set mean age.'], ...
            (1-coverage)*100, n_cpgs_used, n_cpgs_needed, ...
            n_missing_modules + n_singleton_fallback);
    end
end
