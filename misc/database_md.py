"""
If you update the database.csv in any way, this script can be used to update the main
database.md at the top level of this repo.
"""

import pandas as pd

df  = pd.read_csv(f'database.csv')

for col in ['Metabolite name', 'class']:
    df[col] = df[col].str.replace("_", " ", regex=False)

df.fillna("").to_markdown(f'../database.md',
                          index=False,
                          tablefmt='github')
