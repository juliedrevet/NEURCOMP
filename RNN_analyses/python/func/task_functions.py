import numpy as np

def rndvm(n, t, k):
    return np.mod(np.random.vonmises(t * 2, k, n) / 2, np.pi)

def get_seqang(seqtlt, seqcat, catang=None):
    # draw random generative angles corresponding to orange generative category
    if catang is None:
        catang = np.random.uniform(0, np.pi)
    # assign center of blue generative category to second element
    matmap = [catang, catang + np.pi / 2]
    # constrain to the [0,π] interval
    matmap = np.mod(matmap, np.pi)
    # get actual angle corresponding to tilt from generative category
    seqang = np.mod(seqtlt + matmap[seqcat], np.pi)
    return seqang, np.mod(catang - np.pi / 4, np.pi)

def gen_blck(iseqlen, condn=np.arange(3), nseq = 6, kappa = 0.5,sigcat = 10,equalized_per_cond = False):
    # by default: generate 3 conditions (PRE/POST/REF)
    #                    x 2 categories (orange/blue)
    ncond = len(condn)

    if equalized_per_cond:
        seqtlt_base = rndvm([nseq, iseqlen], 0, kappa)
        seqcat_base = np.tile(np.arange(2), nseq // 2 + 1)[:nseq]
        # same tilts and categories for each condition
        seqtlt = np.tile(seqtlt_base, (ncond, 1))  # (nseq*ncond, iseqlen)
        seqcat = np.tile(seqcat_base, ncond)  # (nseq*ncond,)
        seqcon = np.repeat(condn, nseq)  # (nseq*ncond,)
        catang = np.tile(np.random.uniform(0, np.pi, nseq), ncond)
        # noise of the reference (PRE) drawn for all sequences, whatever the conditions: the same seed gives the
        # same sequences when the conditions are simulated separately
        refnoise = np.tile(np.random.normal(0, sigcat * np.pi / 180, size=(nseq, iseqlen)), (ncond, 1))
        total_nseq = nseq * ncond
    else: # no equalization
        seqtlt = rndvm([nseq, iseqlen], 0, kappa)
        seqcat = np.tile(np.arange(2), nseq // 2 + 1)[:nseq]
        seqcon = np.repeat(condn, int(np.ceil(nseq / ncond)))[:nseq]
        catang = None
        total_nseq = nseq

    nseq = total_nseq

    # define input array
    input_array = np.zeros((nseq,iseqlen+1,5))
    # dim2 = iseqlen+1 => last element is post-mapping cue
    # dim3 = 5 => first 2 elements = reference cos and sin , next 2 elements = stimulus cos and sin, last element: GO signal

    seqllr = np.zeros(nseq)
    seqang = [None] * nseq
    refang = np.zeros(nseq)
    cumllr = np.zeros((nseq,iseqlen))
    cumnul = np.zeros((nseq,iseqlen))
    tltllr = np.zeros((nseq,iseqlen))
    tltnul = np.zeros((nseq,iseqlen))

    for iseq in range(nseq):
        cond = seqcon[iseq]

        seq, ref = get_seqang(seqtlt[iseq, :], seqcat[iseq], None if catang is None else catang[iseq])
        seqang[iseq] = seq
        refang[iseq] = ref

        if cond == 0: # PRE
            sigcat_rad = sigcat * np.pi / 180
            if equalized_per_cond:  # drawn above
                ref_noisy = ref + refnoise[iseq]
            else:
                ref_noisy = ref + np.random.normal(0, sigcat_rad, size=iseqlen)

        if seqcat[iseq] == 0:   # generative category = orange category
            tltllr_tmp = 2 * kappa * np.cos(2 * seqtlt[iseq, :])
            tltnul_tmp = 2 * kappa * np.cos(2 * (seqtlt[iseq, :] - np.pi / 4)) # null category
        else: # generative category =  blue category (is orange + pi/2)
            tltllr_tmp = 2 * kappa * np.cos(2 * (seqtlt[iseq, :] + np.pi / 2))
            tltnul_tmp = 2 * kappa * np.cos(2 * (seqtlt[iseq, :] + np.pi / 2 - np.pi / 4))

        seqllr[iseq] = np.sum(tltllr_tmp)
        cumllr[iseq] = np.cumsum(tltllr_tmp)
        cumnul[iseq] = np.cumsum(tltnul_tmp)
        tltllr[iseq] = tltllr_tmp
        tltnul[iseq] = tltnul_tmp

        if cond == 0: # PRE: mapping-cue retrieval from working memory is noisy at each step
            input_array[iseq, :iseqlen, 0] = np.cos(2 * ref_noisy)
            input_array[iseq, :iseqlen, 1] = np.sin(2 * ref_noisy)
        if cond == 2: # REF: mapping-cue is exact
            input_array[iseq, :iseqlen, 0] = np.cos(2 * ref)
            input_array[iseq, :iseqlen, 1] = np.sin(2 * ref)
        if cond == 1 or cond == 2:  # POST and REF: clean reference at last timestep
            input_array[iseq, -1, 0] = np.cos(2 * ref)
            input_array[iseq, -1, 1] = np.sin(2 * ref)

        # stimulus channels
        input_array[iseq, :iseqlen, 2] = np.cos(2 * seq)
        input_array[iseq, :iseqlen, 3] = np.sin(2 * seq)

        # GO signal
        input_array[iseq, -1, -1] = 1

    # define corresponding target array
    category_array = 1 - seqcat

    return (input_array, category_array, seqllr, seqcon, seqang, refang, cumllr, cumnul, tltllr, tltnul)


def generate_batch(batch_size, iseqlen, condn=np.arange(3), sigcat = 10, equalized_per_cond = False):
    ncond = len(condn)

    if equalized_per_cond:
        assert batch_size % ncond == 0, f"batch_size must be divisible by ncond for equalized_per_cond=True"
        nseq = batch_size // ncond  # Output will have exactly batch_size // ncond sequences for each condition
    else:
        nseq = batch_size  # Output will have exactly batch_size sequences

    out = gen_blck(iseqlen, condn, nseq, sigcat=sigcat, equalized_per_cond=equalized_per_cond)

    input_array = out[0]
    category_array = out[1]
    seqllr = out[2]
    seqcon = out[3]
    seqang = out[4]
    refang = out[5]
    cumllr = out[6]
    cumnul = out[7]
    tltllr = out[8]
    tltnul = out[9]

    nseq = len(seqcon)
    if nseq!=batch_size:
        raise ValueError("batch_size not the same as nseq")

    # shuffle generated sequences
    idx_shuffle = np.arange(nseq)
    np.random.shuffle(idx_shuffle)

    seqang_batch = [seqang[i] for i in idx_shuffle]

    # compute optimal performance
    seqcon_tmp = seqcon[idx_shuffle]
    seqllr_tmp = seqllr[idx_shuffle]
    seqcat_tmp = category_array[idx_shuffle]
    perf_opt = np.sum((seqllr_tmp > 0) == (seqcat_tmp == 1)) / len(idx_shuffle)
    perf_opt_cond = np.full(3, np.nan)
    for icond in np.arange(3):
        if np.sum(seqcon_tmp == icond) > 0:
            perf_opt_cond[icond] = np.sum(
                (seqllr_tmp[seqcon_tmp == icond] > 0) == (seqcat_tmp[seqcon_tmp == icond] == 1)) / np.sum(
                seqcon_tmp == icond)

    output_batch = input_array[idx_shuffle, :, :], category_array[idx_shuffle], seqllr[idx_shuffle], seqcon[
        idx_shuffle], perf_opt, perf_opt_cond, cumllr[idx_shuffle, :], cumnul[idx_shuffle, :], tltllr[idx_shuffle, :], tltnul[idx_shuffle, :], refang[idx_shuffle], seqang_batch

    return output_batch