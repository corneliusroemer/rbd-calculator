#%%
import pandas as pd
import numpy as np
from collections import defaultdict
import json

#%%
d = pd.read_csv(
    "https://raw.githubusercontent.com/jbloomlab/SARS2_RBD_Ab_escape_maps/main/processed_data/escape_calculator_data.csv",
    sep=",",
)

#%%
# Get Wuhan Spike
wuhan = """
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFR\
SSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIR\
GWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVY\
SSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQ\
GFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFL\
LKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITN\
LCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCF\
TNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYN\
YLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPY\
RVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFG\
RDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAI\
HADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPR\
RARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTM\
YICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFG\
GFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFN\
GLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQN\
VLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGA\
ISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMS\
ECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAH\
FPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELD\
SFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELG\
KYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSE\
PVLKGVKLHYT\
"""
print(wuhan)

#%%
known_to_neutralize = "Omicron BA.2"
escape_strength = 2
precision = 5

#%%
subset = d.loc[
    d.known_to_neutralize.apply(lambda x: known_to_neutralize in x.split(";"))
].copy()
subset
#%%
# Extract the right data from `;` separated column
subset.apply(
    lambda x: float(
        x.neg_log_IC50.split(";")[
            x.known_to_neutralize.split(";").index(known_to_neutralize)
        ]
    ),
    axis=1,
)
#%%

subset["focal_neg_log_IC50"] = subset.apply(
    lambda x: float(
        x.neg_log_IC50.split(";")[
            x.known_to_neutralize.split(";").index(known_to_neutralize)
        ]
    ),
    axis=1,
)

#%%
# Calculate binding strength remaining in focal?
# Is escape measured in all together and only biased by the focal?
neg_log_IC50 = subset.groupby("condition").focal_neg_log_IC50.mean().to_dict()
max_escape = subset.groupby("condition").escape.max().to_dict()

#%%
escape_data = defaultdict(dict)

for datum in subset.itertuples():
    escape_data[datum.condition][datum.site] = datum.escape

escape_data
#%%
total_weight = sum(list(neg_log_IC50.values()))

nextclade_escape_data = {
    "data": [
        {
            "name": antibody,
            "weight": np.round(
                neg_log_IC50[antibody] / total_weight, precision
            ),
            "locations": {
                # Set to 0 for all Wuhan sites
                str(p - 1): {
                    "default": (
                        -np.round(
                            escape_strength
                            * np.log(max(0.01, 1 - v / max_escape[antibody])),
                            precision,
                        )
                    ),
                    wuhan[p]: 0,
                }
                if p != 493
                else -np.round(
                    escape_strength
                    * np.log(max(0.01, 1 - v / max_escape[antibody])),
                    precision,
                )
                for p, v in by_site.items()
            },
        }
        for antibody, by_site in escape_data.items()
    ]
}

with open("escape_data.json", "w") as fh:
    json.dump(nextclade_escape_data, fh)

# %%
