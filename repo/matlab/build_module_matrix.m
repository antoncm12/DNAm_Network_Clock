function X = build_module_matrix(betas, cpg_ids, W)
%BUILD_MODULE_MATRIX  Aggregate a samples×CpGs matrix into samples×modules.
%
%   X = BUILD_MODULE_MATRIX(betas, cpg_ids, W)
%
%   Inputs
%     betas    [n_samples × n_cpgs]   methylation beta values
%     cpg_ids  {1 × n_cpgs} cell      CpG identifiers matching `betas`
%     W        struct from load_clock_weights
%
%   Output
%     X        [n_samples × n_modules] each column is the per-module value
%              (median of available CpGs for 'median' modules; the single
%              CpG's value for 'singleton' modules; training mean if a
%              module has no input CpGs at all).
%
%   Emits a warning if more than 5% of clock CpGs are absent from `cpg_ids`.

    n_samples = size(betas, 1);
    n_modules = W.n_modules;
    X = zeros(n_samples, n_modules);

    % Index lookup: CpG ID → column in betas (case-sensitive)
    if isstring(cpg_ids), cpg_ids = cellstr(cpg_ids); end
    [~, idx_lookup] = ismember(cpg_ids, cpg_ids);   % identity placeholder
    cpg_to_col = containers.Map(cpg_ids, num2cell(1:numel(cpg_ids)));

    n_missing_modules    = 0;
    n_singleton_fallback = 0;
    n_partial_median     = 0;
    n_cpgs_used          = 0;
    n_cpgs_needed        = 0;

    for j = 1:n_modules
        cpgs     = W.module_to_cpgs{j};
        strategy = W.module_to_strategy{j};
        n_cpgs_needed = n_cpgs_needed + numel(cpgs);

        % Resolve which CpGs are present in the input
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
                if n_present < numel(cpgs)
                    n_partial_median = n_partial_median + 1;
                end
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
             'training-set mean age.'], ...
            (1-coverage)*100, n_cpgs_used, n_cpgs_needed, ...
            n_missing_modules + n_singleton_fallback);
    end
end
