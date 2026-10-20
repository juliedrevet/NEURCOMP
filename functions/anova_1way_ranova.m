function [stat, rm] = anova_1way_ranova(x, name, label, dobf)
%% anova_1way_ranova [claude]
% One-way repeated-measures ANOVA with fitrm/ranova (Statistics and Machine
% Learning Toolbox). The p-value is Greenhouse-Geisser corrected, as in
% repanova_auto (with 2 conditions, epsilon = 1 and nothing changes).
%   x      nsubj x ncond: one row per participant, one column per condition
%   name   name of the factor (optional)                     [default: 'condition']
%   label  title of the printout (optional)
%   dobf   also run a Bayesian repeated-measures ANOVA with bf.anova (bayesFactor
%          toolbox, subject as random effect) (optional)     [default: false]
% Prints and returns:
%   stat.F      F value
%   stat.df     degrees of freedom [effect error] (uncorrected)
%   stat.p      p-value, Greenhouse-Geisser corrected
%   stat.eps    Greenhouse-Geisser epsilon
%   stat.eta2p  partial eta squared
%   stat.bf10   Bayes factor for the effect of the factor (only if dobf = true;
%               printed as BF10 if p < 0.05, as BF01 = 1/bf10 otherwise)
%   rm          repeated-measures model (fitrm), e.g. for post-hoc comparisons:
%               multcompare(rm, name)
%
% Example (replaces repanova_auto(performance,{'condition'}))
%   [stat, rm] = anova_1way_ranova(performance, 'condition');   % performance: nsubj x 3
%   p_corr = multcompare(rm, 'condition');                       % pairwise post-hoc comparisons
if nargin < 2 || isempty(name)
    name = 'condition';
end
if nargin < 3
    label = '';
end
if nargin < 4
    dobf = false;
end
if ~ismatrix(x) || size(x,2) < 2
    error('x must be nsubj x ncond, with at least 2 conditions!');
end
[n, k] = size(x);

% one column per condition, and the within-subject factor
t = array2table(x);
w = table(categorical((1:k)'), 'VariableNames', {name});
rm = fitrm(t, sprintf('%s-%s ~ 1', t.Properties.VariableNames{1}, t.Properties.VariableNames{end}), 'WithinDesign', w);
[tbl, ~, C] = ranova(rm, 'WithinModel', name);

% effect of the factor and its error term
ie = find(strcmp(tbl.Properties.RowNames, ['(Intercept):' name]));
eps_tbl = epsilon(rm, C{2}); % C{1} = intercept
stat.F     = tbl.F(ie);
stat.df    = [tbl.DF(ie) tbl.DF(ie+1)];
stat.p     = tbl.pValueGG(ie);
stat.eps   = eps_tbl.GreenhouseGeisser;
stat.eta2p = tbl.SumSq(ie) / (tbl.SumSq(ie) + tbl.SumSq(ie+1));

% Bayesian repeated-measures ANOVA (as in examples/commonDesigns.m of the
% bayesFactor toolbox); fixed seed because bf.anova uses Monte Carlo integration
if dobf
    lvl = repmat(1:k,n,1);
    sub = repmat((1:n)',1,k);
    tbf = table(x(:),categorical(lvl(:)),categorical(sub(:)),'VariableNames',{'y',name,'subject'});
    s = rng; rng(1);
    stat.bf10 = bf.anova(tbf,sprintf('y~%s+subject',name),'treatAsRandom',{'subject'},'verbose',false) / ...
                bf.anova(tbf,'y~subject','treatAsRandom',{'subject'},'verbose',false);
    rng(s);
end

if ~isempty(label)
    fprintf('%s\n', label);
end
if stat.p < 1e-4
    pstr = '< 0.0001';
else
    pstr = sprintf('= %.4f',stat.p);
end
fprintf('repeated-measures ANOVA (%s, %d levels, n = %d), p Greenhouse-Geisser corrected:\n', name, k, n);
fprintf('\t%s F(%d,%d) = %.4f, p %s (eps = %.4f, partial eta2 = %.4f)\n', name, stat.df(1), stat.df(2), stat.F, pstr, stat.eps, stat.eta2p);
if dobf
    if stat.p < 0.05
        fprintf('\tBayesian repeated-measures ANOVA: BF10 = %.4f\n', stat.bf10);
    else
        fprintf('\tBayesian repeated-measures ANOVA: BF01 = %.4f\n', 1/stat.bf10);
    end
end
end
