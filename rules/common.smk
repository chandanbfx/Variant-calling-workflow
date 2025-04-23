"""
Common functions and wildcard constraints
"""

import pandas as pd

# Load samples at workflow start
try:
    samples_df = pd.read_table(config["sample_sheet"])
    SAMPLES = samples_df["sample"].tolist()
except:
    SAMPLES = list(config.get("samples", {}).keys())

wildcard_constraints:
    sample="|".join(SAMPLES)

# Reference paths (safer access)
try:
    REF_FA = config.get("reference", "reference/genome.fa")
    REF_FAI = f"{REF_FA}.fai"
except NameError:
    REF_FA = "reference/genome.fa"
    REF_FAI = f"{REF_FA}.fai"
