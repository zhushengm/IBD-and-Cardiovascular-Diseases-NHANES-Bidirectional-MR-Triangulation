import sys
sys.path.insert(0, r'C:\Users\zhushengm\AppData\Roaming\Python\Python311\site-packages')

from TwoSampleMR import query
import pandas as pd

print("=== Testing available_outcomes ===")
try:
    df = query.available_outcomes()
    print(f"Total GWAS datasets: {df.shape[0]}")
    print(f"Columns: {list(df.columns)}")
    # Search for perianal abscess
    abscess_df = df[df['trait'].str.lower().str.contains('abscess', na=False)]
    print(f"\nAbscess-related GWAS: {len(abscess_df)}")
    if not abscess_df.empty:
        print(abscess_df[['id','trait','nsnp','sample_size','author','year']].to_string())
    else:
        # Show some example traits
        print("\nSample traits in database:")
        print(df['trait'].head(20).tolist())
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
