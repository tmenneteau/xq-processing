import argparse

import pandas as pd
import numpy as np
from pyteomics import mzxml
from pyteomics import mgf

# This script generates a single MGF file from multiple MGF files using identification information obtained from xQuest

def pep_seq_editing(pep_seq_1, pep_seq_2):
    # Replace Cystein in sequence by Ccm (carbamidomethylation) and X (Oxidated methionine in xQuest) in sequence by Mox for xiVIEW
    pep_seq_1 = pep_seq_1.replace('C', 'Ccm')
    pep_seq_1 = pep_seq_1.replace('X', 'Mox')

    pep_seq_2 = pep_seq_2.replace('C', 'Ccm')
    pep_seq_2 = pep_seq_2.replace('X', 'Mox')
    return pep_seq_1, pep_seq_2

def csv_generator(csv_in, output_filename):
    df = pd.read_csv(csv_in)
    clean_df = []
    columns_df=['PepSeq1','PepSeq2','PepPos1','PepPos2','LinkPos1','LinkPos2','Protein1','Protein2','Charge','CrossLinkerModMass','ScanId','ScanmzXML','PeakListFileName','Original_Peaklist','FragmentTolerance','ExpMz','Score']
    i = 0
    for index, row in df.iterrows():
        # Editing of some parameters to make them xiVIEW friendly
        peptide_sequences = row['id'].split('-')
        pep_seq_1, pep_seq_2 = pep_seq_editing(peptide_sequences[0],peptide_sequences[1])
        pep_pos1 = row['AbsPos1'] - row['top1'] + 1
        pep_pos2 = row['AbsPos2'] - row['top2'] + 1
        precursor_mz = row['mzscans'].split(':')

        # Generation of lists for crosslinked peptides light AND heavy
        light_xl = [pep_seq_1,pep_seq_2,pep_pos1,pep_pos2,row['top1'],row['top2'],row['prot1'],row['prot2'],1,138.06808,i,row['Spectrum_XL_1'],output_filename + '.mgf',row['File'],'20 ppm',precursor_mz[0],row['score']]
        heavy_xl = [pep_seq_1,pep_seq_2,pep_pos1,pep_pos2,row['top1'],row['top2'],row['prot1'],row['prot2'],1,150.14381,i+1,row['Spectrum_XL_2'],output_filename + '.mgf',row['File'],'20 ppm',precursor_mz[1],row['score']]

        # Append to dataframe and incrementing new index for MGF
        clean_df.append(light_xl)
        clean_df.append(heavy_xl)
        i = i + 2
    df_data = pd.DataFrame(clean_df, columns=columns_df)
    df_data_cleaned = df_data.copy()
    df_data_cleaned.drop(['ScanmzXML','Original_Peaklist'], axis=1, inplace=True)
    df_data_cleaned.to_csv(output_filename+'.csv', index=False)
    return df_data

def MGF_generator(df, folder_mzxml, output_filename):
    temp_spectrum_manager = []
    for index, row in df.iterrows():
        fn_mzxml_file = folder_mzxml + row['Original_Peaklist'] + '.mzXML'
        with mzxml.read(fn_mzxml_file) as output :
            spectrum_1 = str(row['ScanmzXML'])
            data_spectrum1 = output[spectrum_1]
            title = str('File:'+row['Original_Peaklist']+ '.' + str(row['ScanmzXML']) + ' "scan=' + str(row['ScanId'])+'"')
            params_dict = {'TITLE': title, 'CHARGE': str('1+'), 'PEPMASS': str(row['ExpMz']), 'SCANS': str(row['ScanId'])}
            dictionnaire = {'params': params_dict, 'm/z array': data_spectrum1['m/z array'], 'intensity array': data_spectrum1['intensity array']}
            temp_spectrum_manager.append(dictionnaire)
            
    output_filename = output_filename + '.mgf'
    mgf.write(temp_spectrum_manager, output_filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generates a unique MGF file from multiple MGF and a csv file output from a modified validate_xl for xiVIEW spectrum visualisation tool')
    parser.add_argument('csv_in', help='CSV file from modified validate_xl')
    parser.add_argument('-mzxml_folder', help='folder with mzxml used for xQuest with a 1-based index', default='./')
    parser.add_argument('-output_filename', help='Name of the peaklist and csv files generated', default='test')
    args = parser.parse_args()

    clean_df = csv_generator(args.csv_in, args.output_filename)
    MGF_generator(clean_df, args.mzxml_folder, args.output_filename)
    print('MGF and CSV generated')