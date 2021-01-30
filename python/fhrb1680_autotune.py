
## Init
import os
import numpy as np
import numpy.random as random
import pandas as pd

from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.csv import CsvDataset
from scvi.inference import UnsupervisedTrainer
from scvi.models import SCANVI, VAE
from scvi.inference.autotune import auto_tune_scvi_model

from umap import UMAP

import torch
import scanpy as sc
import louvain

import logging
import pickle
from hyperopt import hp


# %matplotlib inline

use_cuda = True
n_epochs_all = None
save_path = ''
show_plot = True
os.chdir(path = "/Users/hru/Dropbox/aplastic_anemia/")

## Download data
fhrb1680_1 = CsvDataset(filename='results/scrnaseq/FHRB1680/scvi/fhrb1680_1.csv', save_path='', sep=',', new_n_genes=False)
fhrb1680_2 = CsvDataset(filename='results/scrnaseq/FHRB1680/scvi/fhrb1680_2.csv', save_path='', sep=',', new_n_genes=False)
fhrb1680_3 = CsvDataset(filename='results/scrnaseq/FHRB1680/scvi/fhrb1680_3.csv', save_path='', sep=',', new_n_genes=False)

## Combine
all_dataset = GeneExpressionDataset()
all_dataset.populate_from_per_batch_list(Xs = [fhrb1680_1.X, fhrb1680_2.X, fhrb1680_3.X])

n_epochs        = 50
max_evals       = 100
reserve_timeout = 180
fmin_timeout    = 300

trainer, trials = auto_tune_scvi_model(
    gene_dataset=all_dataset,
    parallel=False,
    use_batches=True,
    exp_key="all_dataset",
    train_func_specific_kwargs={"n_epochs": n_epochs},
    max_evals=max_evals,
    reserve_timeout=reserve_timeout,
    fmin_timeout=fmin_timeout,
)

torch.save(trainer.model.state_dict(), 'results/scrnaseq/FHRB1680/scvi/scvi_best.pkl')

## Sample posterior to get latent representation and save those embeddings
full = trainer.create_posterior(trainer.model, all_dataset, indices=np.arange(len(all_dataset)))

latent, batch_indices, labels = full.sequential().get_latent()
batch_indices = batch_indices.ravel()

np.savetxt("results/scrnaseq/FHRB1680/scvi/batch_latent_best.csv", latent, delimiter=",")
np.savetxt("results/scrnaseq/FHRB1680/scvi/batch_indices_best.csv", batch_indices, delimiter=",")
