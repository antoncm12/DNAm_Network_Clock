function ages = predict_age(betas, cpg_ids, weights_dir)
%PREDICT_AGE  Apply the K=60 DNAm Network Clock to a methylation beta matrix.
%
%   ages = PREDICT_AGE(betas, cpg_ids)
%   ages = PREDICT_AGE(betas, cpg_ids, weights_dir)
%
%   Inputs
%     betas       Numeric matrix of size [n_samples × n_cpgs]. Each ROW is one
%                 sample, each COLUMN is one CpG. NaN values are tolerated.
%     cpg_ids     Cell array of length n_cpgs giving the CpG identifier
%                 (Illumina cg* string) for each column of `betas`. Missing
%                 CpGs are tolerated.
%     weights_dir (Optional) Path to the weights/ folder. If omitted, the
%                 function looks for `../weights/` relative to its own path.
%
%   Output
%     ages        [n_samples × 1] column vector of predicted ages (years),
%                 in the same row order as `betas`.
%
%   Example
%     T = readtable('examples/example_input.csv', 'ReadRowNames', true, ...
%                   'VariableNamingRule', 'preserve');
%     cpgs  = T.Properties.VariableNames;
%     betas = table2array(T);
%     ages  = predict_age(betas, cpgs);
%
%   See also: load_clock_weights, build_module_matrix

% ── Argument handling ────────────────────────────────────────────────────
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
            ['Could not locate weights/ directory automatically. Pass the ' ...
             'path as the third argument.']);
    end
end

if size(betas, 2) ~= numel(cpg_ids)
    error('predict_age:DimMismatch', ...
        'numel(cpg_ids) (%d) must equal number of columns in betas (%d).', ...
        numel(cpg_ids), size(betas, 2));
end

% ── Load weights + build module matrix + apply linear model ─────────────
W = load_clock_weights(weights_dir);
X = build_module_matrix(betas, cpg_ids, W);     % [n_samples × n_modules]

X = bsxfun(@minus,   X, W.scaler_mean');         % standardise
X = bsxfun(@rdivide, X, W.scaler_scale');
X = bsxfun(@minus,   X, W.pca_mean');            % centre for PCA

PC_scores = X * W.pca_components';               % [n_samples × 60]
ages = PC_scores * W.pc_coefs + W.intercept;     % [n_samples × 1]
end
