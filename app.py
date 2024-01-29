import pandas as pd
import matplotlib.pyplot as plt
import PySimpleGUI as sg
import csv

# Global variable to store the DataFrame
df = None

def calculate_analysis(data, reference_gene, control_group_number):
    """
    Calculate the analysis results including ΔCt, ΔΔCt, and fold change for each gene.

    Parameters:
    - data: DataFrame containing the experiment results
    - reference_gene: The internal control gene used for normalization
    - control_group_number: The index of the control group

    Returns:
    - fold_change_dict: Dictionary containing fold change DataFrames for each gene
    """
    # Ensure required columns are present
    required_columns = ['Gene', 'Group', 'Ct']
    missing_columns = set(required_columns) - set(data.columns)
    if missing_columns:
        raise ValueError(f"Required columns {missing_columns} not found in the DataFrame.")

    # Handle missing or non-numeric Ct values
    data['Ct'] = pd.to_numeric(data['Ct'], errors='coerce')

    # Get unique genes excluding the reference gene
    target_genes = data[data['Gene'] != reference_gene]['Gene'].unique()

    fold_change_dict = {}

    for gene in target_genes:
        # Calculate mean Ct for the reference gene in each group
        reference_group = data[data['Gene'] == reference_gene]
        reference_mean = reference_group.groupby(['Group'])['Ct'].mean().reset_index()
        reference_mean.rename(columns={'Ct': f'{reference_gene}_Mean'}, inplace=True)

        # Calculate mean Ct for the target gene in each group
        target_gene_group = data[data['Gene'] == gene]
        target_gene_means = target_gene_group.groupby(['Group'])['Ct'].mean().reset_index()
        target_gene_means.rename(columns={'Ct': f'{gene}_Mean'}, inplace=True)

        # Subtract the reference gene's mean Ct to get ΔCt
        delta_ct = target_gene_means.copy()
        delta_ct[f'{gene}_ΔCt'] = target_gene_means[f'{gene}_Mean'] - reference_mean[f'{reference_gene}_Mean']

        # Calculate ΔΔCt by subtracting control group's ΔCt
        delta_ct[f'{gene}_ΔΔCt'] = delta_ct[f'{gene}_ΔCt'] - delta_ct.loc[delta_ct['Group'] == control_group_number][f'{gene}_ΔCt'].values[0]

        # Calculate fold change using 2^(-ΔΔCt)
        delta_ct[f'{gene}_FoldChange'] = 2 ** -delta_ct[f'{gene}_ΔΔCt']

        # Store the fold change DataFrame in the dictionary
        fold_change_dict[gene] = delta_ct[['Group', f'{gene}_ΔCt', f'{gene}_ΔΔCt', f'{gene}_FoldChange']]

    return fold_change_dict

def write_csv(fold_change_dict):
    with open('result.csv', 'w', newline='', encoding='utf_8_sig') as f:
        writer = csv.writer(f)
        writer.writerow(['Gene', 'Group', 'ΔCt', 'ΔΔCt', 'Fold Change'])
        for gene, df in fold_change_dict.items():
            # Remove 'Gene' column if it exists
            
            if 'Gene' in df.columns:
                df = df.drop(columns=['Gene'])
            df.insert(0, 'Gene', gene)
            writer.writerows(df.values.tolist())

def plot_graphs(fold_change_dict):
    """
    Plot bar graphs for fold change for each gene.

    Parameters:
    - fold_change_dict: Dictionary containing fold change DataFrames for each gene
    """
    for gene, fold_change_df in fold_change_dict.items():
        plt.figure(figsize=(10, 6))
        plt.bar(fold_change_df['Group'], fold_change_df[f'{gene}_FoldChange'], width=0.4, color='black')
        plt.xlabel('Group')
        plt.ylabel('Fold Change')
        plt.title(f'The mRNA expression of {gene}')
        plt.show()

def create_gui():
    sg.theme('LightGrey1')

    layout = [
        [sg.Text("Upload Excel File:")],
        [sg.Input(key='-FILE-', visible=False), sg.FileBrowse(), sg.Button("Read")],
        [sg.Text("Select Reference Gene:"), sg.Combo(values=[], key='-REFERENCE-', readonly=True, size=(15, 1))],
        [sg.Text("Select Control Group:"), sg.Combo(values=[], key='-CONTROL-', readonly=True, size=(15, 1))],
        [sg.Table(values=[], headings=['Gene', 'Group', 'ΔCt', 'ΔΔCt', 'Fold Change'],
                  auto_size_columns=False, justification='right',
                  display_row_numbers=False, num_rows=25, key='-TABLE-')],
        [sg.Button("Analyze"), sg.Button("Plot Graphs"), sg.Exit()]
    ]

    window = sg.Window("CalcDeltaApp", layout, resizable=True)

    global df  # Make df global

    while True:
        event, values = window.read()

        if event == sg.WINDOW_CLOSED or event == "Exit":
            break
        elif event == "Read" and values['-FILE-']:
            try:
                file_path = values['-FILE-']

                # Check if the file is an Excel file
                if not file_path.lower().endswith('.xlsx'):
                    raise ValueError("Please select an Excel file.")

                df = pd.read_excel(file_path, na_values=['N/A'])

                if 'Gene' not in df.columns:
                    raise ValueError("Column 'Gene' not found in the Excel file.")

                gene_list = df['Gene'].unique().tolist()
                window['-REFERENCE-'].update(values=gene_list)
                
                if 'Group' not in df.columns:
                    raise ValueError("Column 'Group' not found in the Excel file.")

                gene_list = df['Group'].unique().tolist()
                window['-CONTROL-'].update(values=gene_list)

            except Exception as e:
                sg.popup_error(f"Error reading Excel file: {e}")

        elif event == "Analyze" and values['-REFERENCE-'] and values['-CONTROL-']:
            try:
                reference_gene = values['-REFERENCE-']
                control_group_number = values['-CONTROL-']

                # Check if the control group number is valid
                if control_group_number not in df['Group'].unique():
                    raise ValueError(f"Invalid control group number: {control_group_number}")

                fold_change_dict = calculate_analysis(df, reference_gene, control_group_number)

                # Display the results in the table
                table_values = []
                for gene, fold_change_df in fold_change_dict.items():
                    # Add gene name to each row
                    fold_change_df.insert(0, 'Gene', gene)
                    table_values.extend(fold_change_df.values.tolist())
                window['-TABLE-'].update(values=table_values)

                # Write results to CSV
                write_csv(fold_change_dict)

            except Exception as e:
                sg.popup_error(f"Error during analysis: {e}")

        elif event == "Plot Graphs" and values['-REFERENCE-'] and values['-CONTROL-']:
            try:
                reference_gene = values['-REFERENCE-']
                control_group_number = values['-CONTROL-']

                # Check if the control group number is valid
                if control_group_number not in df['Group'].unique():
                    raise ValueError(f"Invalid control group number: {control_group_number}")

                fold_change_dict = calculate_analysis(df, reference_gene, control_group_number)
                plot_graphs(fold_change_dict)

            except Exception as e:
                sg.popup_error(f"Error during analysis: {e}")

    window.close()

create_gui()
