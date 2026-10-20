import torch
import numpy as np
from scipy.stats import pearsonr, spearmanr
from scipy.stats import t as tdist
from scipy.integrate import quad

def get_accuracy(logits: torch.Tensor, targets: torch.Tensor):
    pred = torch.round(torch.sigmoid(logits)).detach().cpu().numpy()  # Convert logits to binary predictions
    pcor = np.mean(pred == targets.detach().cpu().numpy())  # Compute accuracy
    return pcor


def get_accuracy_per_cond(logits: torch.Tensor, targets: torch.Tensor, seqcon):
    pred = torch.round(torch.sigmoid(logits)).detach().cpu().numpy()  # Convert logits to binary predictions
    pcor_cond = np.full((3,), float('nan'))
    for icond in range(3):
        if np.sum(seqcon == icond) > 0:
            pcor_cond[icond] = np.mean(pred[seqcon == icond] == targets[seqcon == icond].detach().cpu().numpy())  # Compute accuracy
    return pcor_cond


def ridge_matlab(y, X, k):
    """
    Ridge coefficients as MATLAB's ridge(y,X,k) with default scaling:
    predictors centred and scaled to unit s.d. (N-1), no constant term,
    solution of [Z; sqrt(k)*I] \\ [y; 0] (pseudo-observations).
    y can be (n, nvar): one solution per column of y (same X), b is then (p, nvar).
    """
    n, p = X.shape
    mx = X.mean(axis=0)
    sx = X.std(axis=0, ddof=1)
    sx[np.abs(sx) < np.sqrt(np.finfo(float).eps)] = 1  # as MATLAB: constant predictors not scaled
    Z = (X - mx) / sx
    b = np.linalg.solve(Z.T @ Z + k * np.eye(p), Z.T @ y)
    return b


def decode_ridg(C, D, nfold, alpha, sdtrn=True, corrtype='pearson'):
    """
    Decode variable using ridge regression (as decode_ridg.m).

    Parameters
    ----------
    C        : (n,) variable to decode
               or (n, nvar): several variables decoded from the same data D
               (same folds, one decoder per variable), outputs then have one column per variable
    D        : (n, k) data matrix (e.g. hidden units)
    nfold    : number of cross-validation folds (interleaved: observation j, 1-based,
               belongs to test fold mod(j, nfold)), or list of (train, test) index arrays
    alpha    : ridge regularization factor
    sdtrn    : fold-wise standardization of the data w.r.t. the training set
    corrtype : 'pearson' or 'spearman'

    Returns
    -------
    dict with Ctru (ground truth), Chat (decoded value), mhat, shat (fold-wise data
    mean and s.d.), What (fold-wise weights), rdec (decoding correlation coefficient
    on all test predictions pooled), Wdec (overall weights without cross-validation)
    for C (n, nvar): Ctru, Chat (n, nvar), What (nfeat, nvar, nfold), rdec (nvar,), Wdec (nfeat, nvar)
    """
    C = np.asarray(C, dtype=float)
    is1d = C.ndim == 1
    C = C.reshape(C.shape[0], -1)
    D = np.asarray(D, dtype=float)
    if D.ndim != 2 or D.shape[0] != C.shape[0]:
        raise ValueError('Mismatching input sizes!')
    n, nfeat = D.shape
    nvar = C.shape[1]

    # use or create training and test sets
    if isinstance(nfold, (list, tuple)):
        itrn = [np.asarray(f[0]).ravel() for f in nfold]
        itst = [np.asarray(f[1]).ravel() for f in nfold]
        nfold = len(nfold)
    else:
        j = np.arange(1, n + 1)  # 1-based observation number, as in MATLAB
        itst = [np.flatnonzero(j % nfold == ifold) for ifold in range(nfold)]
        itrn = [np.flatnonzero(j % nfold != ifold) for ifold in range(nfold)]

    # get decoding set
    idec = np.concatenate(itst)
    if np.unique(idec).size < idec.size:
        raise ValueError('Repeated observations found in decoding set!')

    # perform decoding analysis
    Chat = np.full((n, nvar), np.nan)
    mhat = np.full((nfold, nfeat), np.nan)
    shat = np.full((nfold, nfeat), np.nan)
    What = np.full((nfeat, nvar, nfold), np.nan)
    for ifold in range(nfold):
        Dhat = D
        if sdtrn:
            # standardize data w.r.t. training set
            mhat[ifold] = D[itrn[ifold]].mean(axis=0)
            shat[ifold] = D[itrn[ifold]].std(axis=0, ddof=1)
            Dhat = (D - mhat[ifold]) / shat[ifold]
        # compute decoding weights on training set
        W = ridge_matlab(C[itrn[ifold]], Dhat[itrn[ifold]], alpha)
        # compute predictions on test set
        Chat[itst[ifold]] = Dhat[itst[ifold]] @ W
        What[:, :, ifold] = W

    if corrtype == 'pearson':
        rdec = np.array([pearsonr(C[idec, ivar], Chat[idec, ivar])[0] for ivar in range(nvar)])
    elif corrtype == 'spearman':
        rdec = np.array([spearmanr(C[idec, ivar], Chat[idec, ivar])[0] for ivar in range(nvar)])
    else:
        raise ValueError(f'unknown correlation type {corrtype}')

    # compute overall decoding weights on standardized data
    Dz = (D - D.mean(axis=0)) / D.std(axis=0, ddof=1)
    Wdec = ridge_matlab(C[idec], Dz[idec], alpha)

    if is1d:
        C, Chat, What, rdec, Wdec = C[:, 0], Chat[:, 0], What[:, 0, :], rdec[0], Wdec[:, 0]
    return {'Ctru': C, 'Chat': Chat, 'mhat': mhat, 'shat': shat,
            'What': What, 'rdec': rdec, 'Wdec': Wdec}


def cluster_permutation(Y, X=None, alpha_form=0.05, tail='both', nperm=10000, minsize=2, seed=1234,
                        xvec=None, verbose=True, Z=None):
    """
    Cluster-based permutation test on time courses (Maris & Oostenveld, 2007): adjacent
    samples whose t-value exceeds the cluster-forming threshold are grouped into clusters,
    and the tstat of each cluster (sum of its t-values) is compared to the null distribution
    of the maximum absolute cluster tstat across permutations.

    Parameters
    ----------
    Y          : (n, nsmp) one row per participant/network (e.g. time course, or epoch course)
    X          : (n,) variable to correlate with Y (Pearson, permutation of the rows), or None
                 to test Y against 0 (one-sample t-test, random sign flips of the rows).
                 For a paired test, pass the difference between the two conditions as Y.
                 X can also be (n, nsmp): one variable per sample, correlated with the
                 same sample of Y (the rows are permuted the same way for all samples).
    alpha_form : cluster-forming threshold (uncorrected p of each sample)
    tail       : 'both', 'right' (positive effects) or 'left'
    nperm      : number of permutations
    minsize    : minimum number of adjacent samples in a cluster
    seed       : seed of the random number generator (own generator, global state unchanged)
    xvec       : (nsmp,) value of each sample (time, epoch...), for printing
    verbose    : print the clusters
    Z          : (n,) or (n, nsmp) covariate for a partial correlation between X and Y:
                 the linear effect of Z (with intercept) is removed from X and Y at each sample,
                 the residuals of X are permuted, df = n - 3

    Returns
    -------
    dict with t (nsmp,), r (nsmp,, correlation only), df, tcrit, clusters (list of dict:
    idx, tstat, p, xbeg, xend) and null_tstat (nperm,), with the cluster-corrected
    p = (1 + #(null_tstat >= |tstat|)) / (nperm + 1)
    """
    Y = np.asarray(Y, dtype=float)
    n, nsmp = Y.shape
    if xvec is None:
        xvec = np.arange(1, nsmp + 1)
    rng = np.random.default_rng(seed)

    if X is not None:  # correlation with X across rows
        X = np.asarray(X, dtype=float).reshape(n, -1)
        if X.shape[1] not in (1, nsmp):
            raise ValueError('X must have one value per row of Y (or one column per sample of Y)!')
        df = n - 2
        if Z is not None:
            Z = np.asarray(Z, dtype=float).reshape(n, -1)
            X, Y = _residuals(X, Z), _residuals(Y, Z)
            df = n - 3
        if X.shape[1] == 1:  # same X for all samples
            X = X[:, 0]
            r = _corr(X[:, None], Y)[0]
            # permutations: shuffle rows
            Xperm = np.stack([X[rng.permutation(n)] for _ in range(nperm)], axis=1)  # (n, nperm)
            rperm = _corr(Xperm, Y)                                                   # (nperm, nsmp)
        else:
            r = _corr_paired(X, Y)
            rperm = np.stack([_corr_paired(X[rng.permutation(n)], Y) for _ in range(nperm)])  # (nperm, nsmp)
        t = r * np.sqrt(df / (1 - r**2))
        tperm = rperm * np.sqrt(df / (1 - rperm**2))
    else:  # one-sample t-test against 0
        df = n - 1
        t = Y.mean(axis=0) / (Y.std(axis=0, ddof=1) / np.sqrt(n))
        # permutations: flip the sign of rows
        S = 2 * (rng.random((n, nperm)) > 0.5) - 1                 # (n, nperm)
        m = S.T @ Y / n                                            # (nperm, nsmp)
        v = (np.sum(Y**2, axis=0) - n * m**2) / (n - 1)            # squares do not depend on the signs
        tperm = m / np.sqrt(v / n)

    # cluster-forming threshold
    if tail == 'both':
        tcrit, sgn = tdist.ppf(1 - alpha_form / 2, df), (+1, -1)
    elif tail == 'right':
        tcrit, sgn = tdist.ppf(1 - alpha_form, df), (+1,)
    elif tail == 'left':
        tcrit, sgn = tdist.ppf(1 - alpha_form, df), (-1,)
    else:
        raise ValueError(f"unknown tail '{tail}': use 'both', 'right' or 'left'")

    # observed clusters and null distribution of the maximum absolute cluster tstat
    cidx, ctstat = _find_clusters(t, tcrit, sgn, minsize)
    null_tstat = np.zeros(nperm)
    for iperm in range(nperm):
        _, tstat_perm = _find_clusters(tperm[iperm], tcrit, sgn, minsize)
        if tstat_perm:
            null_tstat[iperm] = np.max(np.abs(tstat_perm))

    clusters = []
    for idx, ts in zip(cidx, ctstat):
        p = (1 + np.sum(null_tstat >= abs(ts))) / (nperm + 1)
        clusters.append({'idx': idx, 'tstat': ts, 'p': p, 'xbeg': xvec[idx[0]], 'xend': xvec[idx[-1]]})

    if verbose:
        print(f'{len(clusters)} cluster(s) found')
        for c in clusters:
            print(f"\tcluster tstat (sum) : {c['tstat']:+5.3f} on [{c['xbeg']} {c['xend']}]")
        for ic, c in enumerate(clusters):
            print(f'Comparing cluster {ic + 1}/{len(clusters)} :')
            print(f"\tcluster-corrected p = {c['p']:.5f}")

    out = {'t': t, 'df': df, 'tcrit': tcrit, 'clusters': clusters, 'null_tstat': null_tstat}
    if X is not None:
        out['r'] = r
    return out


def _corr(A, Y):
    # Pearson correlation between each column of A (n, na) and each column of Y (n, nsmp): (na, nsmp)
    A = (A - A.mean(axis=0)) / A.std(axis=0)
    Y = (Y - Y.mean(axis=0)) / Y.std(axis=0)
    return A.T @ Y / A.shape[0]


def _corr_paired(X, Y):
    # Pearson correlation between each column of X (n, nsmp) and the same column of Y: (nsmp,)
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    Y = (Y - Y.mean(axis=0)) / Y.std(axis=0)
    return np.mean(X * Y, axis=0)


def _residuals(A, Z):
    # residuals of the linear regression (with intercept) of each column of A on the
    # matching column of Z (A and Z: (n, 1) or (n, nsmp))
    A, Z = np.broadcast_arrays(A, Z)
    A = A - A.mean(axis=0)
    Z = Z - Z.mean(axis=0)
    return A - Z * np.sum(Z * A, axis=0) / np.sum(Z**2, axis=0)


def _find_clusters(t, tcrit, sgn, minsize):
    # clusters of at least minsize adjacent samples with s*t > tcrit, for each sign s in sgn
    cidx, ctstat = [], []
    for s in sgn:
        above = np.concatenate([[False], s * t > tcrit, [False]])
        ibeg = np.flatnonzero(np.diff(above.astype(int)) == 1)
        iend = np.flatnonzero(np.diff(above.astype(int)) == -1) - 1
        for b, e in zip(ibeg, iend):
            if e - b + 1 >= minsize:
                cidx.append(np.arange(b, e + 1))
                ctstat.append(np.sum(t[b:e + 1]))
    return cidx, ctstat


def bf_ttest(t, n, scale=np.sqrt(2) / 2):
    """
    BF10 of a one-sample (or paired) t-test from its t-value and the number of observations n (df = n - 1)
    """
    df = n - 1

    def integrand(g):
        return ((1 + n * g * scale**2) ** -0.5 * (1 + t**2 / ((1 + n * g * scale**2) * df)) ** (-(df + 1) / 2)
                * (2 * np.pi) ** -0.5 * g ** -1.5 * np.exp(-1 / (2 * g)))

    num = quad(integrand, 0, np.inf)[0]
    den = (1 + t**2 / df) ** (-(df + 1) / 2)
    return num / den


def partial_corr(X, Y, Z):
    # partial correlation between X and Y controlling for Z, column-wise (rows = networks):
    # correlation between the residuals of the linear regressions (with intercept) of X and Y on Z
    X, Y, Z = np.broadcast_arrays(*(np.asarray(a, dtype=float) for a in (X, Y, Z)))
    Z = Z - Z.mean(axis=0)
    Xr = X - X.mean(axis=0)
    Yr = Y - Y.mean(axis=0)
    Xr = Xr - Z * np.sum(Z * Xr, axis=0) / np.sum(Z**2, axis=0)
    Yr = Yr - Z * np.sum(Z * Yr, axis=0) / np.sum(Z**2, axis=0)
    return np.sum(Xr * Yr, axis=0) / np.sqrt(np.sum(Xr**2, axis=0) * np.sum(Yr**2, axis=0))