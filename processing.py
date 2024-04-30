import pandas as pd

# Step 1: Read the CSV file containing bioactivity data
df = pd.read_csv('CALU1_raw.csv')

# Step 2: Filter compounds with bioactivity data (IC50, EC50, or GI50)
df = df[df['Standard Type'].isin(['IC50', 'EC50', 'GI50'])]

# Step 3: Unit conversion to μM
def convert_to_uM(row):
    if row['Standard Units'] == 'nM':
        return row['Standard Value'] / 1000  # Convert nM to μM
    elif row['Standard Units'] == 'M':
        return row['Standard Value'] * 1000  # Convert M to μM
    elif row['Standard Units'] == 'ug/mL':
        # Additional conversion logic if needed
        pass
    else:
        return row['Standard Value']  # Assume already in μM

df['standard_value_uM'] = df.apply(convert_to_uM, axis=1)

# Step 4: Calculate average bioactivity value for molecules with multiple records
df_avg = df.groupby('Molecule ChEMBL ID')['standard_value_uM'].mean().reset_index()

# Step 5: Assign activity labels based on bioactivity values
def assign_activity_label(value):
    if value <= 10:
        return 'active'
    else:
        return 'inactive'

df_avg['activity_label'] = df_avg['standard_value_uM'].apply(assign_activity_label)

# Step 6: Exclude molecules with ambiguous labels
df_filtered = df_avg[((df_avg['activity_label'] == 'active') | ((df_avg['activity_label'] == 'inactive') & (df_avg['standard_value_uM'] >= 50)))]


# Merge the filtered DataFrame with the original DataFrame to include 'Smiles' column
df_processed = pd.merge(df_filtered, df[['Molecule ChEMBL ID', 'Smiles']], on='Molecule ChEMBL ID', how='left')

# Save the processed data to a new CSV file
df_processed.to_csv('processed_data_with_smiles.csv', index=False)
