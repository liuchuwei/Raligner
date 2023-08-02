# Load VEGA and packages
import sys
import vega
import pandas as pd
import scanpy as sc
# Ignore warnings in this tutorial
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import torch
import random
import os

def setup_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True

setup_seed(888)

if not os.path.exists("tcga/tmp"):
    # Loading and preparing data
    tcga_ccle = sc.read_loom("data/tcga_cell.loom")
    # adata = vega.data.pbmc()
    print(tcga_ccle)

    vega.utils.setup_anndata(tcga_ccle)

    # Creating and training the VEGA model
    model = vega.VEGA(tcga_ccle,
                      gmt_paths='reactomes.gmt',
                      add_nodes=1,
                      positive_decoder=True)

    print(model)

    # save model
    model.save('./tcga/tmp', save_adata=True, save_history=True)


for id in range(98,100):
    # out path
    out_path = "tcga/latent_"+ str(id) + ".csv"

    # load model
    model = vega.VEGA.load('./tcga/tmp')
    model.train_vega(n_epochs=50, use_gpu=True)

    # output latent
    latent = model.to_latent()
    df = pd.DataFrame(latent)
    df.columns = model.adata.uns['_vega']['gmv_names']
    df.to_csv(out_path, index=False)

    # gene weight
    w = model.decoder._get_weights().data
    w = pd.DataFrame(np.array(w.cpu()))

    # out path
    out_path = "tcga/weight_"+ str(id) + ".csv"

    w.to_csv(out_path, index=False)


