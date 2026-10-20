import torch
import torch.nn as nn

class CustomRNN(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, noise_value=0.0, xavier_gain= 1.0):
        super().__init__()

        self.hidden_size = hidden_size
        self.noise_value = noise_value
        self.sel_noise_dir = None
        self.sel_noise_value = 0.0
        self.output_size = output_size
        self.xavier_gain = xavier_gain

        self.i2h = nn.Linear(input_size, hidden_size, bias=False)  # input to hidden
        self.h2h = nn.Linear(hidden_size, hidden_size, bias=True)  # hidden to hidden
        self.h2o = nn.Linear(hidden_size, output_size, bias=False)  # hidden to output

        # standard init
        gain = self.xavier_gain * nn.init.calculate_gain("tanh")

        nn.init.xavier_normal_(self.i2h.weight, gain=gain)
        nn.init.orthogonal_(self.h2h.weight)
        nn.init.xavier_normal_(self.h2o.weight)

        nn.init.zeros_(self.h2h.bias)

        self.tanh = nn.Tanh()

    def forward(self, input, hidden, device, add_noise=True, rec_rnd=None, sel_rnd=None):
        # rec_rnd (batch, hidden_size) and sel_rnd (batch, 1): standard-normal samples of the recurrent and of
        # the selective noise, given by test (the same samples for a sequence in each condition), drawn here if None

        inp = self.i2h(input)
        hid = self.h2h(hidden)

        combined = inp + hid

        if add_noise and (self.noise_value != 0):
            if rec_rnd is None:
                noise_samp = torch.normal(mean=torch.zeros(combined.shape, device=device),
                                          std=torch.ones(combined.shape, device=device) * self.noise_value)
            else:
                noise_samp = rec_rnd * self.noise_value
            combined += noise_samp

        if add_noise and (self.sel_noise_value != 0):  # selective noise along sel_noise_dir: one sample per sequence and stimulus
            if sel_rnd is None:
                sel_rnd = torch.randn(combined.shape[0], 1, device=device)
            combined = combined + sel_rnd * self.sel_noise_value * self.sel_noise_dir

        hidden = self.tanh(combined)

        logits = self.h2o(hidden).squeeze(-1)
        return logits, hidden, combined

    def initHidden(self, batch_size, device):
        # initialize hidden state
        return torch.zeros(batch_size, self.hidden_size, device=device)

