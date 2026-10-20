import os
import torch.optim as optim
from func.CustomRNN import *
from func.task_functions import *
from func.get_metrics import *


def train(model, batch_size, nbatch, condn_train, device, optimizer, seq2train=range(2,13), scheduler=None, sigcat = 10):
    # One epoch loops over nbatch gradient updates.
    # Each gradient update sees ALL sequence lengths (stacked before loss).
    # Loss is computed against the generative category (target).

    model.train()
    nseqlen = len(seq2train)

    for batch_idx in range(nbatch):

        perf_opt = np.zeros(nseqlen)
        perf_opt_cond = np.zeros((nseqlen, 3))
        seqcon = np.zeros((batch_size, nseqlen))

        # accumulate outputs and targets across ALL sequence lengths before loss
        all_outputs = []
        all_targets = []
        all_taropt  = []

        for iseqlen, jseqlen in enumerate(seq2train):
            # generate sequences of length jseqlen+1
            out_train = generate_batch(batch_size, jseqlen, condn_train, sigcat = sigcat)
            data = out_train[0]

            perf_opt[iseqlen]       = out_train[4]
            perf_opt_cond[iseqlen]  = out_train[5]
            seqcon[:, iseqlen]      = out_train[3]

            # targets: generative category label
            target_sl = torch.tensor(out_train[1], dtype=torch.float32, device=device)
            # taropt: optimal choice given observed sequence (ceiling = 100%)
            taropt_sl = torch.tensor(np.sign(out_train[2]), dtype=torch.float32, device=device).clamp(min=0)

            data   = torch.from_numpy(data).float().to(device)
            ncues  = data.size(1)  # seqlen + 1 (stimuli + POST-cue/GO)

            hidden = model.initHidden(batch_size, device)
            for i in range(ncues):
                if i < ncues - 1:
                    add_noise = True
                elif i == ncues - 1:
                    add_noise = False
                logits, hidden, _ = model(data[:, i, :], hidden, device,add_noise=add_noise)
            # logits after last timestep (POST-cue + GO signal)

            all_outputs.append(logits)
            all_targets.append(target_sl)
            all_taropt.append(taropt_sl)

        # stack across sequence lengths: (batch_size, nseqlen)
        output = torch.stack(all_outputs, dim=1)
        target = torch.stack(all_targets, dim=1)
        taropt = torch.stack(all_taropt,  dim=1)

        # single gradient update seeing all sequence lengths
        loss = torch.nn.functional.binary_cross_entropy_with_logits(output, target)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        if scheduler:
            scheduler.step()

    # metrics from last batch only
    pcor      = get_accuracy(output, target)
    popt      = get_accuracy(output, taropt)
    pcor_cond = get_accuracy_per_cond(output, target, seqcon)
    popt_cond = get_accuracy_per_cond(output, taropt, seqcon)

    return pcor, pcor_cond, loss.item(), np.mean(perf_opt), np.mean(perf_opt_cond,axis=0), popt, popt_cond


def run_training(irnn,noise_value,cfg,verbose = False):
    # early stopper if patience > 0
    # default adam optimizer

    input_size = cfg['input_size']
    hidden_size = cfg['hidden_size']
    output_size = cfg['output_size']
    alpha = cfg['alpha']
    batch_size = cfg['batch_size']
    condn_train = cfg['condn_train']
    nbatch = cfg['nbatch']
    max_epochs = cfg['max_epochs']
    fname_train = cfg['fname_train']
    device = cfg['device']
    saveornot = cfg['saveornot']
    patience = cfg['patience']
    # set default cfg values
    if 'sigcat' in cfg:
        sigcat = cfg['sigcat']
    else:
        sigcat = 10 # in ° degrees not rad
    if 'xavier_gain' in cfg:
        xavier_gain = cfg['xavier_gain']
    else:
        xavier_gain = 1.0
    if 'seqlen' in cfg:
        seq_lens = cfg['seqlen']
    else:
        seq_lens = range(2,13)  # sequence lengths: 2 to 12
    if 'savecheck' in cfg:
        savecheck = cfg['savecheck']
    else:
        savecheck = False

    seed = int(np.random.SeedSequence(cfg.get('seed_base', 0), spawn_key=(irnn,)).generate_state(1)[0])
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.set_num_threads(1)

    fname_noise = fname_train + (f'_noise_{noise_value:.2f}'.replace('.', '-') if noise_value != 0 else '')
    os.makedirs(fname_noise, exist_ok=True)
    fname_irnn = os.path.join(fname_noise, f'model_irnn_{irnn:02d}')

    model = CustomRNN(input_size, hidden_size, output_size, noise_value, xavier_gain = xavier_gain)

    model.to(device)
    optimizer = optim.Adam(model.parameters(), lr=alpha)

    history_train = {'loss': [], 'acc': [], 'opt': [], 'accopt': [],
                     'PRE': [], 'POST': [], 'REF': [],
                     'PREopt': [], 'POSTopt': [], 'REFopt': [],
                     'PREaccopt': [], 'POSTaccopt': [], 'REFaccopt': []}

    if patience >0:
        early_stopper = EarlyStoppingByAccuracy(threshold=0.95, patience=patience)

    if savecheck>0:
        torch.save(model.state_dict(), f'{fname_irnn}_epo_0000.pt')

    for epoch in range(max_epochs):

        out_train = train(model, batch_size, nbatch, condn_train, device, optimizer, seq_lens, sigcat = sigcat)

        accuracy = out_train[0] * 100
        accuracy_cond = out_train[1] * 100
        loss = out_train[2]
        perf_opt = out_train[3] * 100
        perf_opt_cond = out_train[4] * 100
        accuracy_opt = out_train[5] * 100
        accuracy_opt_cond = out_train[6] * 100

        history_train['loss'].append(loss)
        history_train['acc'].append(accuracy)
        history_train['opt'].append(perf_opt)
        history_train['accopt'].append(accuracy_opt)
        history_train['PRE'].append(accuracy_cond[0])
        history_train['POST'].append(accuracy_cond[1])
        history_train['REF'].append(accuracy_cond[2])
        history_train['PREopt'].append(perf_opt_cond[0])
        history_train['POSTopt'].append(perf_opt_cond[1])
        history_train['REFopt'].append(perf_opt_cond[2])
        history_train['PREaccopt'].append(accuracy_opt_cond[0])
        history_train['POSTaccopt'].append(accuracy_opt_cond[1])
        history_train['REFaccopt'].append(accuracy_opt_cond[2])

        if verbose:
            if (max_epochs <= 1000) or (max_epochs >= 1000 and epoch % 20 == 0):
                print(f'[Epoch {epoch + 1}/{max_epochs} ] loss: {loss:.4f}, acc: {accuracy:2.2f}%, accopt: {accuracy_opt:2.2f}%')

        if patience >0 and early_stopper.check(accuracy_opt):
            print(f"Early stopping triggered at epoch {epoch + 1}")
            break

        if savecheck>0 and (epoch+1) % savecheck == 0:
            torch.save(model.state_dict(), f'{fname_irnn}_epo_{epoch+1:04d}.pt')

    rnn_output = {'max_epochs': float(max_epochs),
                  'condition_train': condn_train.astype(float),
                  'nbatch': float(nbatch),
                  'batch_size': float(batch_size),
                  'learning_rate': alpha,
                  'noise_value': float(noise_value),
                  'seqlen': seq_lens,
                  'seed': float(seed),
                  'history_train': history_train}

    for name, param in model.named_parameters():
        rnn_output[name.replace('.', '_')] = param.detach().numpy()

    if saveornot:
        print(f'trained RNN saved as {fname_irnn}')
        torch.save({
                 'epoch': epoch,
                 'model_state_dict': model.state_dict(),
                 'optimizer_state_dict': optimizer.state_dict(),
                 'loss': loss,
                 'acc': accuracy,
                 'cfg': cfg,
                 **rnn_output
             }, f'{fname_irnn}.tar')

    return history_train, model, optimizer

class EarlyStoppingByAccuracy:
    def __init__(self, threshold=0.95, patience=400):
        self.threshold = threshold
        self.patience = patience
        self.counter = 0

    def check(self, acc):
        if acc >= self.threshold*100:
            self.counter += 1
        else:
            self.counter = 0
        return self.counter >= self.patience

    def reset(self):
        self.counter = 0