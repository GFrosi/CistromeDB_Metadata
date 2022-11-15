import pandas as pd


def filter_human(df):
   '''Receives a df (Cistrome metadata).
    Returns a df filtered by Human and 
    ChIP-Seq data''' 
   
   return df[(df['species'].str.contains('Homo sapiens', case=False, na=False))\
         & (~df['target'].str.contains('DNAse|-SEQ', na=False, case=False))]



def filter_human_noENC(df_human):
    '''Receives a df (Human ChIP-Seq
    Cistrome metadata). Returns a df
    filtered by ENCODE samples''' 

    return df_human[~df_human['GSM'].str.contains('ENC', na=False)]


