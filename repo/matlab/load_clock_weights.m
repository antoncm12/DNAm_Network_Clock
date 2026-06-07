function W = load_clock_weights(weights_dir)
%LOAD_CLOCK_WEIGHTS  Load the K=60 DNAm Network Clock weights from CSV files.
%
%   W = LOAD_CLOCK_WEIGHTS(weights_dir)
%
%   Returns a struct W with the following fields, all in canonical module
%   column order (length 4595 unless otherwise stated):
%     .module_order        [4595×1 double]  module IDs
%     .module_to_cpgs      {4595×1 cell}    CpG IDs per module
%     .module_to_strategy  {4595×1 cell}    'median' or 'singleton'
%     .scaler_mean         [4595×1 double]
%     .scaler_scale        [4595×1 double]
%     .pca_mean            [4595×1 double]
%     .pca_components      [60×4595 double] PC loadings (rows = PCs)
%     .pc_coefs            [60×1 double]    ridge coefficients
%     .intercept           scalar           ridge intercept
%     .K                   scalar = 60
%     .n_modules           scalar = 4595
%
%   Results are cached in a persistent variable so repeated calls don't
%   re-read the CSVs.

    persistent cache cache_key
    key = weights_dir;
    if ~isempty(cache) && strcmp(key, cache_key)
        W = cache; return;
    end

    % ── module_definitions.csv ──────────────────────────────────────────
    defs_path = fullfile(weights_dir, 'module_definitions.csv');
    if ~isfile(defs_path)
        error('load_clock_weights:NoFile', 'Missing %s', defs_path);
    end
    defs = readtable(defs_path, 'TextType', 'string');

    % ── module_params.csv ──────────────────────────────────────────────
    params = readtable(fullfile(weights_dir, 'module_params.csv'));
    module_order = params.Module;
    n_modules    = numel(module_order);

    % Per-module CpG lists, in module_order
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

    % ── pca_components_k60.csv ─────────────────────────────────────────
    comps_tbl = readtable(fullfile(weights_dir, 'pca_components_k60.csv'), ...
                          'VariableNamingRule', 'preserve');
    % Drop the 'PC' label column → numeric matrix [60 × 4595]
    pca_components = table2array(comps_tbl(:, 2:end));
    if size(pca_components, 1) ~= 60 || size(pca_components, 2) ~= n_modules
        error('load_clock_weights:BadComponents', ...
            'pca_components_k60.csv has shape %dx%d, expected 60x%d.', ...
            size(pca_components, 1), size(pca_components, 2), n_modules);
    end

    % ── ridge_k60.csv ──────────────────────────────────────────────────
    ridge = readtable(fullfile(weights_dir, 'ridge_k60.csv'), 'TextType', 'string');
    pc_mask   = startsWith(ridge.name, "PC");
    pc_coefs  = ridge.value(pc_mask);
    intercept = ridge.value(ridge.name == "intercept");
    if numel(intercept) ~= 1
        error('load_clock_weights:BadRidge', ...
            'ridge_k60.csv must contain exactly one row with name="intercept".');
    end

    % ── Assemble and cache ────────────────────────────────────────────
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
