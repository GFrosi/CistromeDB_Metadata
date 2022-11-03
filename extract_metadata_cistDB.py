from cgitb import reset
import pandas as pd
import sys
import requests
import pandas as pd
from tqdm import tqdm
from time import sleep
import logging
from retry import retry
import logging
from utils.loggerinitializer import *
from distutils.dir_util import mkpath
from json import JSONDecodeError
import json
import os
import argparse
from bs4 import BeautifulSoup


mkpath(os.getcwd() + "/logs/")
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
initialize_logger(os.getcwd() + "/logs/", logger)


@retry(TimeoutError, tries=5, delay=3)
def dict_of_dict(number, path):
    """
    Receives a number to create a range
    to build a url. Returns a dict of dict
    where each key is the number used to build
    the url and the dict values contains the 
    the metadata information
    """

    logger.info("You are accessing " + str(number) + "pages")
    #number = int(number)
    dict_master = {}
    root = "http://dc2.cistrome.org/api/inspector?id="
    myfile = open(os.path.join(path,'data'+str(number)+'.json'), 'w')

    logger.info("Starting requests...")
    for num in tqdm(range(number, 105615)):
        
        dict_partial = {}
        dict_partial[num] = {}
        
        try:
            dict_master[num] = {} #creating key 
            url = root+str(num) #buinding url
            
            #dict_data = requests.get(url).json() #inbuilt json() constructor, dict object
            response = requests.get(url)
            sleep(2)
            data = response.text
            soup = BeautifulSoup(data, 'lxml').find_all('body')
            for i in soup:
                if "Request denied" not in i.get_text():          
                    dict_partial[num] = {}
                    dict_data = response.json()
                    # print(dict_data)       
                    
                    #getting metadata info
                    for ele in dict_data["treats"]: #ele is a dict
                        dict_master[num]['GSE'] = 'GSE'+ele['other_ids'].split(":")[-1].strip().replace('"','').replace('}', '') #ele[other_ids] is a string. e.g {"sra": "18393", "gse": "200017937"}
                    
                    dict_master[num]["cell_line"] = dict_data["treats"][0].get('cell_line__name','----') #the list 'treats' has one element (a dict)
                    dict_master[num]["target"] = dict_data["treats"][0].get('factor__name','----')
                    dict_master[num]["is_correcting"] = dict_data["treats"][0].get('is_correcting','----')
                    dict_master[num]["paper_reference"] = dict_data["treats"][0].get('paper__reference','----')
                    dict_master[num]["paper_pmid"] = dict_data["treats"][0].get('paper__pmid','----')
                    dict_master[num]["paper_journal_name"] = dict_data["treats"][0].get('paper__journal__name','----')
                    dict_master[num]["name"] = dict_data["treats"][0].get('name','----')
                    dict_master[num]["disease_state_name"] = dict_data["treats"][0].get('disease_state__name','----')
                    dict_master[num]["disease_state_name"] = dict_data["treats"][0].get('disease_state__name','----')
                    dict_master[num]["cell_type"] = dict_data["treats"][0].get('cell_type__name','----')
                    dict_master[num]["GSM"] = dict_data["treats"][0].get('unique_id','----')
                    dict_master[num]["species"] = dict_data["treats"][0].get('species__name','----')
                    dict_master[num]["tissue_type"] = dict_data["treats"][0].get('tissue_type__name','----')
                    dict_master[num]["paper_lab"] = dict_data["treats"][0].get('paper__lab','----')
        
                    #getting sign
                    dict_master[num]["sign"] = dict_data.get('sign', '----')

                    #Getting QC field
                    #getting QC info judge
                    dict_master[num]["judge_map"] = dict_data.get("qc",{}).get("judge",{}).get('map', '----')
                    dict_master[num]["judge_peaks"] = dict_data.get("qc",{}).get("judge",{}).get('peaks', '----')
                    dict_master[num]["judge_fastqc"] = dict_data.get("qc",{}).get("judge",{}).get('fastqc', '----')
                    dict_master[num]["judge_frip"] = dict_data.get("qc",{}).get("judge",{}).get('frip', '----') 
                    dict_master[num]["judge_pbc"] = dict_data.get("qc",{}).get("judge",{}).get('pbc', '----')
                    dict_master[num]["judge_motif"] = dict_data.get("qc",{}).get("judge",{}).get('motif_judge', '----')
                    dict_master[num]["judge_dhs"] = dict_data.get("qc",{}).get("judge",{}).get('dhs', '----')

                    #getting QC info table
                    # dict_master[num]["map_perc"] = ' '.join(map(str, dict_data["qc"]["table"]["map"]))
                    dict_master[num]["map_perc"] = ','.join(map(str, dict_data.get('qc',{}).get('table',{}).get('map','----')))
                    dict_master[num]["peaks"] = ','.join(map(str, dict_data.get('qc',{}).get('table',{}).get('peaks','----')))
                    dict_master[num]["control_number"] =  dict_data.get("qc",{}).get("table",{}).get("control_number", '----')
                    dict_master[num]["fastqc"] = ''.join(map(str, dict_data.get('qc',{}).get('table',{}).get('fastqc','----')))
                    dict_master[num]["frip"] = ''.join(map(str, dict_data.get('qc',{}).get('table',{}).get('frip','----'))) #getting list results as strings
                    dict_master[num]["sample"] =  ''.join(map(str, dict_data.get('qc',{}).get('table',{}).get('sample','----')))
                    dict_master[num]["meta"] = ''.join(map(str,dict_data.get('qc',{}).get('table',{}).get('meta','----')))
                    dict_master[num]["map_number"] = ''.join(map(str, dict_data.get('qc',{}).get('table',{}).get('map_number','----')))
                    dict_master[num]["pbc_qc"] = ''.join(map(str, dict_data.get('qc',{}).get('table',{}).get('pbc','----')))
                    dict_master[num]["dhs"] = ''.join(map(str, dict_data.get('qc',{}).get('table',{}).get('dhs','----')))
                    dict_master[num]["raw"] = ''.join(map(str, dict_data.get('qc',{}).get('table',{}).get('raw_number','----')))
                    
                    #dict_partial
                    dict_partial[num] = dict_master[num]
                    dict_sep = dict() #separator to be replaced (}{}{ by ,)
                    json.dump(dict_partial, myfile)
                    json.dump(dict_sep, myfile)

                
        except JSONDecodeError as jd:
            logger.error("JSONDecodeError. Restart the program from: " + str(num) + " Something wrong with the url.")
            print(f'JSONDecodeError {jd}. Restart the program from {num}. Something wrong with the url')
            sys.exit(1)

        except requests.exceptions.Timeout:
            logger.error("Timeout! You should restart the program from: " + str(num))
            print(f'Timeout! You should restart the program from: {num}')
            sys.exit(1)

        except ConnectionAbortedError as cae:
            logger.error("Connection aborted error! You shold restart the program from: " + str(num))
            print(f"Error {cae}. You shold restart the program  from {num}")
            sys.exit(1)
                
        except ConnectionRefusedError as cre:
            logger.error("Connection aborted error! You shold restart the program from: " + str(num))
            print(f"Error {cre}. You shold restart the program  from {num}")
            sys.exit(1)
            
    logger.info("Resquests successfully finished!")
    print("Resquests successfully finished!")
    
    return dict_master


def combine_json(path_out): #to work 
    """Receives a path with json files generated
    by the function above. Returns a big json file
    generated from the other json files"""

    list_files = [fname for fname in os.listdir(path_out)]
    result = []

    for f1 in list_files:
        with open(os.path.join(path_out,f1), 'r') as infile:
            result.extend(json.load(infile))

    with open('complete.json', 'w') as output_file:
        json.dump(result, output_file)


def create_df(dict_master): #working
    """Receives a dict of dict
    and returns a df where the
    keys of the first dict are
    the index ang the keys of
    the dict values are the rows"""
    
    df = pd.DataFrame.from_dict(dict_master,orient='index')
    
    return df

  
def main():

    logger.info("Starting program! You have " + str(args.num) + " requests to do!")
    print(f"Starting program! You have {args.num} requests to do!")
    dict_master = dict_of_dict(args.num, args.path)
    logger.info("Requestes finished. Creating dataframe...")
    print("Requestes finished. Creating dataframe...")
    df = create_df(dict_master) #beautiful!
    df.to_csv(args.output)
    logger.info("Dataframe saved!")
    print("Dataframe saved!")


    # if args.concat:
    #     logger.info("Running concat jsons!")
    #     dict_master = combine_json(args.path)
    #     print(dict_master)
    
    #     # create_df(dict_master) #beautiful!
    #     sys.exit()
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description = 'A script to get the the metadata information.'
    )

    parser.add_argument('-n', '--num', action='store',
                        help='The number to create a range to build the urls',
                        type=int, required=True)
    
    parser.add_argument('-c', '--concat',
                        help='The absolute path to the txt files generated in the previous run to be concatenated and mapped into a final df',
                        type=bool, default=False, required=False)

    parser.add_argument('-p', '--path', action='store',
                        help='The absolute path to save the txt file if needed',
                        required=True)
    parser.add_argument('-o', '--output', action='store',
                        help='The absolute path to save the final df',
                        required=True)


    args = parser.parse_args()

    main()
