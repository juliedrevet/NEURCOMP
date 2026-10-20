# simulations of the networks
#   test               test of a network on batches of sequences: same steps as train (training.py), without
#                      gradient (eval mode, torch.no_grad), returning hidden activity and responses after each stimulus
#   sim_rnn            one simulation of a network (one epoch, one noise setting): hidden activity, accuracy
#   dec_rnn            decoding of the variables from the hidden activity of a simulation
#   run_simulations    simulations and decoding for several settings: epochs, recurrent noise, input noise,
#                      noise along a direction (one simulation per value, a single value applies to all simulations)
#   get_directions     decoding direction of each variable (decoder weights), for the causal perturbation
#   run_perturbation   causal perturbation: noise along the decoding direction of each variable
import numpy as np
import torch
from scipy.stats import pearsonr
from func.CustomRNN import CustomRNN
from func.task_functions import *
from func.get_metrics import *


def test(model, batch_size, nbatch, condn_test, device, seq2test=range(2, 13), sigcat=10, var2keep=(), seed=None,
         input_noise=0):
    # test of a network on nbatch batches of batch_size sequences per number of stimuli of seq2test: same steps as
    # train, without gradient. Each sequence is tested in each condition of condn_test, with the same noise samples
    # (recurrent, selective and input noise)
    #   var2keep    : variables returned after each stimulus (cumllr, cumnul, cumsumcos, cumsumsin, tltllr)
    #   seed        : seed of the sequences and of the noise samples (always drawn: the same seed gives the same
    #                 samples whatever the noise levels)
    #   input_noise : s.d. (degrees) of the noise on the orientation of each stimulus, input only (the variables use
    #                 the noise-free stimuli)
    # returns one dict per number of stimuli of seq2test, all batches pooled (N = nbatch x batch_size sequences):
    #   hidden   (N, hidden_size, seqlen) hidden activity after each stimulus
    #   combined (N, hidden_size, seqlen) pre-activation (with noise) after each stimulus
    #   var      dict: value of each variable of var2keep after each stimulus, (N, seqlen)
    #   response (N,) response of the network: 1 = category 1, 0 = category 2
    #   target   (N,) generative category
    #   seqllr   (N,) log-likelihood ratio of the sequence (category 1 vs 2)
    #   seqcon   (N,) condition
    if seed is not None:
        np.random.seed(seed)
        torch.manual_seed(seed)

    ncond = len(condn_test)
    if batch_size % ncond != 0:
        raise ValueError('batch_size must be a multiple of the number of conditions')
    nseq = batch_size // ncond  # sequences per condition

    model.eval()
    out = [{'hidden': [], 'combined': [], 'var': {v: [] for v in var2keep}, 'response': [], 'target': [], 'seqllr': [],
            'seqcon': []} for _ in seq2test]
    with torch.no_grad():
        for batch_idx in range(nbatch):
            for iseqlen, seqlen in enumerate(seq2test):
                # the same nseq sequences in each condition (one condition after the other)
                out_test = gen_blck(seqlen, np.asarray(condn_test), nseq, sigcat=sigcat, equalized_per_cond=True)
                var = {'cumllr':    out_test[6],
                       'cumnul':    out_test[7],
                       'tltllr':    out_test[8],
                       'cumsumcos': np.cumsum(out_test[0][:, :-1, 2], axis=1),
                       'cumsumsin': np.cumsum(out_test[0][:, :-1, 3], axis=1)}

                # noise samples (standard normal) of each sequence, the same in each condition
                rec_rnd = torch.randn(nseq, seqlen, model.hidden_size, device=device).repeat(ncond, 1, 1)
                sel_rnd = torch.randn(nseq, seqlen, 1, device=device).repeat(ncond, 1, 1)
                inp_rnd = torch.randn(nseq, seqlen, device=device).repeat(ncond, 1)

                data = torch.from_numpy(out_test[0]).float().to(device)
                if input_noise != 0:  # noise on the orientation of each stimulus (input only)
                    ang2 = torch.atan2(data[:, :-1, 3], data[:, :-1, 2])  # 2 x orientation
                    ang2 = ang2 + 2 * inp_rnd * input_noise * np.pi / 180
                    data[:, :-1, 2] = torch.cos(ang2)
                    data[:, :-1, 3] = torch.sin(ang2)

                hidden = model.initHidden(batch_size, device)
                hidden_all = torch.zeros((batch_size, model.hidden_size, seqlen), device=device)
                combined_all = torch.zeros((batch_size, model.hidden_size, seqlen), device=device)
                for i in range(seqlen):  # stimuli
                    _, hidden, combined = model(data[:, i, :], hidden, device, rec_rnd=rec_rnd[:, i],
                                                sel_rnd=sel_rnd[:, i])
                    hidden_all[:, :, i] = hidden
                    combined_all[:, :, i] = combined
                logits, _, _ = model(data[:, -1, :], hidden, device, add_noise=False)  # POST-cue, no noise

                o = out[iseqlen]
                o['hidden'].append(hidden_all.cpu().numpy())
                o['combined'].append(combined_all.cpu().numpy())
                for v in var2keep:
                    o['var'][v].append(var[v])
                o['response'].append(torch.round(torch.sigmoid(logits)).cpu().numpy())
                o['target'].append(out_test[1])
                o['seqllr'].append(out_test[2])
                o['seqcon'].append(out_test[3])

    # all batches pooled
    for o in out:
        for k in ('hidden', 'combined', 'response', 'target', 'seqllr', 'seqcon'):
            o[k] = np.concatenate(o[k])
        o['var'] = {v: np.concatenate(x) for v, x in o['var'].items()}
    return out


def sim_rnn(file_epo, cfg, seq2test, nseq, seed, device, condn_test=(0, 1), noise_value=0,
            var2keep=('cumllr', 'cumnul', 'cumsumcos', 'cumsumsin'), sel_noise_dir=None, sel_noise_value=0,
            input_noise=0):
    # test the network saved at one epoch (file_epo, .pt) on equalized sequences (each sequence is
    # tested in each condition) and pool all stimuli of all sequences
    #   seq2test    : numbers of stimuli of the test sequences (e.g. [4, 8, 12] as for humans)
    #   nseq        : number of sequences per number of stimuli (each one tested in each condition)
    #   seed        : seed of the test sequences and of the noise: the same seed gives the same sequences, e.g. for
    #                 all epochs of a network
    #   noise_value : s.d. of the recurrent noise (0: exact network)
    #   sel_noise_dir, sel_noise_value : direction (hidden_size,) and s.d. of a selective noise (see CustomRNN)
    #   input_noise : s.d. (degrees) of the noise on the orientation of each stimulus, input only (see test)
    # returns one dict per condition of condn_test:
    #   hidden     (nobs, hidden_size) hidden activity after each stimulus of each sequence
    #   combined   (nobs, hidden_size) pre-activation (with noise) after each stimulus of each sequence
    #   var        dict: value of each variable of var2keep after each stimulus, (nobs,)
    #   correct    (nobs,) the response to the sequence is its generative category
    #   pcor       accuracy (fraction of sequences)
    #   pcor_bayes accuracy of the Bayes-optimal observer on the same sequences: the sign of the LLR of the sequence
    #              gives its generative category
    model = CustomRNN(cfg['input_size'], cfg['hidden_size'], cfg['output_size'], noise_value)
    model.load_state_dict(torch.load(file_epo, weights_only=True))
    model.to(device)
    if sel_noise_dir is not None:  # selective noise along a direction
        model.sel_noise_dir = torch.as_tensor(sel_noise_dir, dtype=torch.float32, device=device)
        model.sel_noise_value = sel_noise_value

    out_test = test(model, nseq * len(condn_test), 1, np.asarray(condn_test), device, seq2test, cfg['sigcat'],
                    var2keep, seed, input_noise)
    sim = []
    for cond in condn_test:
        s = {'hidden': [], 'combined': [], 'var': {v: [] for v in var2keep}, 'correct': [], 'corseq': [], 'bayseq': []}
        for seqlen, out in zip(seq2test, out_test):
            m = out['seqcon'] == cond
            corseq = out['response'][m] == out['target'][m]  # correct response
            # (nseq, hidden, seqlen) -> (nseq*seqlen, hidden): the stimuli of each sequence one after the other
            s['hidden'].append(out['hidden'][m].transpose(0, 2, 1).reshape(-1, cfg['hidden_size']))
            s['combined'].append(out['combined'][m].transpose(0, 2, 1).reshape(-1, cfg['hidden_size']))
            for v in var2keep:
                s['var'][v].append(out['var'][v][m].reshape(-1))
            s['correct'].append(np.repeat(corseq, seqlen))
            s['corseq'].append(corseq)
            s['bayseq'].append((out['seqllr'][m] > 0) == (out['target'][m] == 1))  # correct response of the optimal observer
        sim.append({'hidden':     np.concatenate(s['hidden']),
                    'combined':   np.concatenate(s['combined']),
                    'var':        {v: np.concatenate(x) for v, x in s['var'].items()},
                    'correct':    np.concatenate(s['correct']),
                    'pcor':       np.mean(np.concatenate(s['corseq'])),
                    'pcor_bayes': np.mean(np.concatenate(s['bayseq']))})
    return sim


def dec_rnn(sim, var2decode, nfold=10, alpha=500):
    # decode each variable from the hidden activity of all stimuli pooled (decode_ridg, as for
    # the MEG data), separately in each condition, and on the stimuli of the correct and of the error
    # sequences (same cross-validated predictions)
    #   sim : output of sim_rnn
    # returns rdec, rdec_cor, rdec_err: dict var -> (ncond,) decoding correlation coefficient (r)
    ncond = len(sim)
    rdec, rdec_cor, rdec_err = ({v: np.full(ncond, np.nan) for v in var2decode} for _ in range(3))
    for icond, s in enumerate(sim):
        cor = s['correct']
        # all variables decoded at once (same hidden activity and folds, one decoder per variable)
        dec = decode_ridg(np.column_stack([s['var'][v] for v in var2decode]), s['hidden'], nfold, alpha)
        for ivar, v in enumerate(var2decode):
            Ctru, Chat = dec['Ctru'][:, ivar], dec['Chat'][:, ivar]
            rdec[v][icond] = dec['rdec'][ivar]
            if np.sum(cor) >= 3:
                rdec_cor[v][icond] = pearsonr(Ctru[cor], Chat[cor])[0]
            if np.sum(~cor) >= 3:
                rdec_err[v][icond] = pearsonr(Ctru[~cor], Chat[~cor])[0]
    return rdec, rdec_cor, rdec_err


def run_simulations(files_epo, cfg, seq2test, nseq, seed, device, var2decode, condn_test=(0, 1), noise_values=0,
                    input_noises=0, sel_noise_dir=None, sel_noise_values=0, nfold=10, alpha=500):
    # simulate a network (sim_rnn) and decode (dec_rnn) for each setting, all on the same sequences (same seed):
    #   files_epo        : file(s) of the network saved at given epochs (.pt), e.g. all epochs (along training)
    #   noise_values     : s.d. of the recurrent noise (with the same seed, the same noise samples are only scaled)
    #   input_noises     : s.d. (degrees) of the noise on the orientation of each stimulus
    #   sel_noise_dir    : direction (hidden_size,) of a selective noise, the same for all simulations (None: no
    #                      selective noise)
    #   sel_noise_values : s.d. of the selective noise
    # one simulation per element of the arguments given as lists or arrays (all of the same length); an argument
    # given as a single value applies to all simulations
    # returns dict with pcor (nsim, ncond), pcor_bayes (accuracy of the Bayes-optimal observer on the same sequences),
    # rdec, rdec_cor, rdec_err: dict var -> (nsim, ncond)
    torch.set_num_threads(1)  # one thread per network (networks simulated in parallel)
    settings = {'file_epo': files_epo, 'noise_value': noise_values, 'input_noise': input_noises,
                'sel_noise_value': sel_noise_values}
    settings = {k: [v] if isinstance(v, str) or np.ndim(v) == 0 else list(v) for k, v in settings.items()}
    nsim = max(len(v) for v in settings.values())
    for k, v in settings.items():
        if len(v) == 1:
            settings[k] = v * nsim
        elif len(v) != nsim:
            raise ValueError(f'{k}: {len(v)} values for {nsim} simulations')
    ncond = len(condn_test)
    out = {'pcor': np.zeros((nsim, ncond)),
           'rdec': {v: np.zeros((nsim, ncond)) for v in var2decode},
           'rdec_cor': {v: np.zeros((nsim, ncond)) for v in var2decode},
           'rdec_err': {v: np.zeros((nsim, ncond)) for v in var2decode}}
    for isim in range(nsim):
        sim = sim_rnn(settings['file_epo'][isim], cfg, seq2test, nseq, seed, device, condn_test,
                      settings['noise_value'][isim], var2decode, sel_noise_dir=sel_noise_dir,
                      sel_noise_value=settings['sel_noise_value'][isim], input_noise=settings['input_noise'][isim])
        out['pcor'][isim] = [s['pcor'] for s in sim]
        dec = dec_rnn(sim, var2decode, nfold, alpha)
        for key, d in zip(('rdec', 'rdec_cor', 'rdec_err'), dec):
            for v in var2decode:
                out[key][v][isim] = d[v]
    out['pcor_bayes'] = sim[0]['pcor_bayes']  # same sequences in all conditions and for all simulations
    return out


def get_directions(file_epo, cfg, seq2test, nseq, seed, device, var2perturb, cond, nfold=10, alpha=500):
    # decoding direction of each variable of var2perturb, in condition cond, from a
    # simulation without noise: decoding weights of all sequences (Wdec of decode_ridg)
    # returns w (hidden_size, nvar), each column normalized
    sim = sim_rnn(file_epo, cfg, seq2test, nseq, seed, device, [cond], 0, var2perturb)[0]
    dec = decode_ridg(np.column_stack([sim['var'][p] for p in var2perturb]), sim['hidden'], nfold, alpha)
    w = dec['Wdec']
    return w / np.linalg.norm(w, axis=0)


def run_perturbation(file_epo, cfg, sel_noise_array, seq2test, nseq, seed_dir, seed_test, device, var2decode,
                     var2perturb=('cumllr', 'cumnul', 'cumsumcos', 'cumsumsin'), condn_test=(0, 1), nfold=10, alpha=500):
    # causal perturbation: noise injected along the decoding direction of each
    # variable (decoder weights), separately in each condition
    #   1. on a first set of sequences (seed_dir), no noise: decoding direction w of each variable of var2perturb
    #      (get_directions)
    #   2. on a second set of sequences (seed_test), no noise: accuracy and decoding (baseline), and s.d. of
    #      the pre-activation along w (scale)
    #   3. on the second set, for each variable of var2perturb and each noise level s: noise s * scale * z * w
    #      added to the pre-activation after each stimulus (z ~ N(0,1), one sample per sequence and stimulus,
    #      the same samples for all variables and levels), accuracy and decoding (decoders retrained, run_simulations)
    # returns dict with pcor0 (ncond,), rdec0 {var: (ncond,)} (baseline), pcor {pvar: (nlev, ncond)},
    # rdec {pvar: {var: (nlev, ncond)}}, w {pvar: (ncond, hidden_size)}, scale {pvar: (ncond,)}
    torch.set_num_threads(1)  # one thread per network (networks tested in parallel)
    ncond, nlev = len(condn_test), len(sel_noise_array)
    out = {'pcor0': np.zeros(ncond), 'rdec0': {v: np.zeros(ncond) for v in var2decode},
           'pcor': {p: np.zeros((nlev, ncond)) for p in var2perturb},
           'rdec': {p: {v: np.zeros((nlev, ncond)) for v in var2decode} for p in var2perturb},
           'w': {p: np.zeros((ncond, cfg['hidden_size'])) for p in var2perturb},
           'scale': {p: np.zeros(ncond) for p in var2perturb}}
    for icond, cond in enumerate(condn_test):
        # 1. directions (first set of sequences)
        w = get_directions(file_epo, cfg, seq2test, nseq, seed_dir, device, var2perturb, cond, nfold, alpha)
        # 2. baseline (second set of sequences)
        sim = sim_rnn(file_epo, cfg, seq2test, nseq, seed_test, device, [cond], 0, var2decode)
        out['pcor0'][icond] = sim[0]['pcor']
        rdec = dec_rnn(sim, var2decode, nfold, alpha)[0]
        for v in var2decode:
            out['rdec0'][v][icond] = rdec[v][0]
        # 3. perturbations (second set of sequences)
        for ip, p in enumerate(var2perturb):
            scale = np.std(sim[0]['combined'] @ w[:, ip], ddof=1)
            out['w'][p][icond] = w[:, ip]
            out['scale'][p][icond] = scale
            res = run_simulations(file_epo, cfg, seq2test, nseq, seed_test, device, var2decode, [cond],
                                  sel_noise_dir=w[:, ip], sel_noise_values=np.asarray(sel_noise_array) * scale,
                                  nfold=nfold, alpha=alpha)
            out['pcor'][p][:, icond] = res['pcor'][:, 0]
            for v in var2decode:
                out['rdec'][p][v][:, icond] = res['rdec'][v][:, 0]
    return out
