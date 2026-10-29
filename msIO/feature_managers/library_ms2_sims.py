import logging

import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

from msIO.feature_managers.db import Library
from msIO.metrics import cosine_similarity_sym


min_sim = .3
max_dmz_da_ms2 = 10e-3
lib = Library(r"C:\Users\Yannick Zander\Downloads\library_high_confidence.sqlite")

logging.basicConfig(level=logging.INFO)
ms2s = lib.get_ms_spectra(lib.feature_ids, level=2)

n_features = len(lib.feature_ids)

sims = []

for f_id1 in tqdm(lib.feature_ids, desc='calculating tanimoto similarities', total=len(lib.feature_ids), smoothing=1/50):
    for f_id2 in lib.feature_ids:
        if f_id2 >= f_id1:
            continue
        if (f_id2 not in ms2s) or (f_id1 not in ms2s):
            continue
        if (scr := cosine_similarity_sym(ms2s[f_id1], ms2s[f_id2], max_dmz_da_ms2)[0]) > min_sim:
            sims.append((f_id1, f_id2, scr))

sims = pd.DataFrame(sims, columns=['f_id_1', 'f_id_2', 'tanimoto_sim'])

G = nx.from_pandas_edgelist(sims, source='f_id_1', target='f_id_2')

nx.draw(G, with_labels=True)
plt.show()

