function [stat, rm] = anova_2way_ranova(x, names, label)
%% anova_2way_ranova [claude]
% Two-way repeated-measures ANOVA with fitrm/ranova (Statistics and Machine
% Learning Toolbox). p-values are Greenhouse-Geisser corrected, as in
% repanova_auto (for factors with 2 levels, epsilon = 1 and nothing changes).
%   x      nsubj x n1 x n2: dimension 2 = factor names{1} (n1 levels),
%          dimension 3 = factor names{2} (n2 levels)
%   names  names of the two factors, in the order of the dimensions of x,
%          e.g. {'position','soundtype'}
%   label  title of the printout (optional)
% Prints and returns stat(e).effect, .F, .df (uncorrected), .p (Greenhouse-Geisser),
% .eps (Greenhouse-Geisser epsilon) and .eta2p (partial eta squared) for
% e = 1: main effect of names{1}, 2: main effect of names{2}, 3: interaction,
% and the repeated-measures model rm (fitrm), e.g. for post-hoc comparisons:
%   multcompare(rm, names{1})                  % levels of names{1}, averaged over names{2}
%   multcompare(rm, names{1}, 'By', names{2})  % levels of names{1}, within each level of names{2}
%
% Example (replaces repanova_auto(sigma,{'seqlen','condition'}))
%   [stat, rm] = anova_2way_ranova(sigma, {'seqlen','condition'});   % sigma: nsubj x 3 x 3
if nargin < 3
    label = '';
end
if ndims(x) ~= 3
    error('x must be nsubj x n1 x n2!');
end
[n, n1, n2] = size(x);

% one column per cell, names{1} varying fastest (as reshape)
t = array2table(reshape(x, n, n1*n2));
w = table(categorical(repmat((1:n1)',n2,1)), categorical(kron((1:n2)',ones(n1,1))), 'VariableNames', names);
rm = fitrm(t, sprintf('%s-%s ~ 1', t.Properties.VariableNames{1}, t.Properties.VariableNames{end}), 'WithinDesign', w);
[tbl, ~, C] = ranova(rm, 'WithinModel', sprintf('%s*%s', names{:}));

effect = {names{1}, names{2}, sprintf('%s x %s',names{1},names{2})};
rows = {sprintf('(Intercept):%s',names{1}), sprintf('(Intercept):%s',names{2}), sprintf('(Intercept):%s:%s',names{1},names{2})};
if ~isempty(label)
    fprintf('%s\n', label);
end
fprintf('repeated-measures ANOVA (%s x %s, n = %d), p Greenhouse-Geisser corrected:\n', names{1}, names{2}, n);
stat = struct('effect',effect,'F',NaN,'df',NaN,'p',NaN,'eps',NaN,'eta2p',NaN);
for e = 1:3
    ie = find(strcmp(tbl.Properties.RowNames, rows{e}));
    stat(e).F   = tbl.F(ie);
    stat(e).df  = [tbl.DF(ie) tbl.DF(ie+1)]; % effect and its error term
    stat(e).p   = tbl.pValueGG(ie);
    eps_tbl     = epsilon(rm, C{e+1});        % C{1} = intercept
    stat(e).eps = eps_tbl.GreenhouseGeisser;
    stat(e).eta2p = tbl.SumSq(ie) / (tbl.SumSq(ie) + tbl.SumSq(ie+1)); % partial eta squared
    if stat(e).p < 1e-4
        pstr = '< 0.0001';
    else
        pstr = sprintf('= %.4f',stat(e).p); 
    end
    fprintf('\t%-*s F(%d,%d) = %.4f, p %s (eps = %.4f, partial eta2 = %.4f)\n', max(cellfun(@numel,effect)), effect{e}, ... 
        stat(e).df(1), stat(e).df(2), stat(e).F, pstr, stat(e).eps, stat(e).eta2p);
end
end
