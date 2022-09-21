#%%
import pandas as pd

#%%
df = pd.read_csv("https://media.githubusercontent.com/media/jbloomlab/SARS-CoV-2-RBD_DMS_Omicron/main/results/final_variant_scores/final_variant_scores.csv")
df
# %%
ba2 = df[df.target == "Omicron_BA2"]
# %%
ba2.set_index(["position", "mutant"], inplace=True)
ba2

#%%
wuhan = df[df.target == "Wuhan-Hu-1_v1"]
wuhan.set_index(["position", "mutant"], inplace=True)

# %%

#%%
def score_mutation(df, mutations):
    df = df.copy()
    rbd_delta = 0
    for pos, aa in mutations:
        # TODO: ignore reversions to Wuhan
        if aa == wuhan.loc[(pos, "A")].wildtype:
            continue
        rbd_delta += df.loc[(pos, aa), "delta_bind"]
    return rbd_delta
#%%
# Example mutations for which we want to calculate the impact on binding
mutations = [(331, "N"), (346, "T"), (460, "K"), (501, "Y")]
score_mutation(ba2, mutations)
# %%
