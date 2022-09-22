#%%
import pandas as pd
import json

#%%
df = pd.read_csv("https://media.githubusercontent.com/media/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron/main/results/final_variant_scores/final_variant_scores.csv")
df
# %%
ba2 = df[df.target == "Omicron_BA2"]
ba2.set_index(["position", "mutant"], inplace=True)
ba2
#%%
ba1 = df[df.target == "Omicron_BA1"]
ba1.set_index(["position", "mutant"], inplace=True)
ba1
#%%
wuhan = df[df.target == "Wuhan-Hu-1_v1"]
wuhan.set_index(["position", "mutant"], inplace=True)
#%%
scores = {}
for pos1, aa in ba2.index:
    pos0 = pos1 - 1
    if pos0 not in scores:
        scores[pos0] = {}
    # if aa == wuhan.loc[(pos, "A")].wildtype:
    #     delta = 0.0
    # else:
    delta = ba2.loc[(pos1, aa), "delta_bind"]
    if pd.isnull(delta):
        delta = ba1.loc[(pos1, aa), "delta_bind"]
    scores[pos0][aa] = delta
scores
#%%
json.dump(scores, open("rbd_scores.json", "w"), indent=2)
# %%
