
import os
import json
import sys

workflow_path = os.environ.get('path_workflow')
if(workflow_path==None):
    raise Exception('You forgot to setup the path_workflow environment variable with the path to the workflow')
    
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

sys.path.insert(0, workflow_path)

from data_processing import ProcessPathwayScores
from build_scoring_matrix import BuildScoringMatrix
from drug_ranking_analysis import DrugRankingAnalysis

class Pipeline_drugResponseCalibration:
    
    def _validate_input(self, config):
        flag = False
        if( os.path.isfile(config) ):
            try:
                with open(config, 'r') as g:
                    cfg = json.load(g)
                    
                aux = True
                if('experiments' in cfg):
                    for e in cfg['experiments']:
                        if( ('folder' in e) and ('expression_file' in e) and ('identifier' in e) ):
                            if ( (e['folder']!='' and e['folder']!=None) and (e['expression_file']!='' and e['expression_file']!=None) and (e['identifier']!='' and e['identifier']!=None) ):
                                if( os.path.isdir(e['folder']) ):
                                    folder = e['folder']
                                    ide = e["identifier"]
                                    if(folder[-1]=='/'):
                                        folder = folder[:-1]
                                    if( os.path.isfile(folder+'/'+e['expression_file']) ):
                                        aux = aux and True
                                    else:
                                        aux = False
                                        print (f'Error - {ide}: Expression file was not found in {folder}')
                                else:
                                    aux = False
                                    print(f'Error - {ide}: The directory was not found in {folder}')
                            else:
                                aux = False
                                print("Error - Mandatory fields are empty")
                        else:
                            aux = False
                            print("Error - Configuration file does not have the mandatory fields")
                else:
                    aux = False
                    print("Error - Configuration file in a bad format")
            except:
                aux = False
                print("Error - Configuration file cannot not be loaded and parsed")
        else:
            print("Error - Configuration file not found")
        
        flag = aux
           
        return flag
    
    def _validate_labels(self, e):
        flag = False
        if( 'labels_file' in e ):
            if ( (e['labels_file']!='' and e['labels_file']!=None) ):
                folder = e['folder']
                ide = e["identifier"]
                if(folder[-1]=='/'):
                    folder = folder[:-1]
                if( os.path.isfile(folder+'/'+e['labels_file']) ):
                    flag = True
                else:
                    print (f'Error - {ide}: Expression file was not found in {folder}')
            else:
                print("Error - Label File field is empty")
        else:
            print("Error - Label File field is missing from configuration")
            
        return flag
    
    def _check_norm(self, e):
        norm = 'fpkm_uq'
        gtf = ''
        
        flag = False
        if( ('type_normalization' in e) and ('genome_assembly' in e) ):
            if ( (e['type_normalization']!='' and e['type_normalization']!=None) and (e['genome_assembly']!='' and e['genome_assembly']!=None) ):
                if( (e['type_normalization'] in ['tpm', 'tmm', 'cpm', 'fpkm', 'fpkm_uq']) and (e['genome_assembly'] in ['g37', 'g38']) ):
                    norm = e['type_normalization']
                    gtf = e['genome_assembly']
                else:
                    if( not (e['type_normalization'] in ['tpm', 'tmm', 'cpm', 'fpkm', 'fpkm_uq']) ):
                        print("Information - type_normalization is not among available options (tpm, fpkm, tmm, cpm, fpkm_uq), reverting to fpkm_uq")    
                    if( not (e['genome_assembly'] in ['g37', 'g38']) ):
                        print("Information - genome_assembly is not among available options (g37, g38), reverting to fpkm_uq")    
            else:
                print("Information - type_normalization and/or genome_assembly are empty, reverting to fpkm_uq")
        else:
            print("Information - type_normalization and/or genome_assembly are missing, reverting to fpkm_uq ")
            
        return norm, gtf
    
    def _check_geneset(self, e):
        path = 'KEGG_2021_HUMAN'
        pok = open( f'{workflow_path}/genesets_available.txt','r').read().split('\n')
        
        flag = False
        if( ('pathway_geneset' in e) ):
            if ( (e['pathway_geneset']!='' and e['pathway_geneset']!=None) ):
                if( (e['pathway_geneset'] in pok) ):
                    path = e['pathway_geneset']
                else:
                    print("Information - pathway_geneset is not among available options in genesets_available.txt, reverting to KEGG_2021_HUMAN")    
            else:
                print("Information - pathway_geneset is empty, reverting to KEGG_2021_HUMAN")
        else:
            print("Information - pathway_geneset is missing, reverting to KEGG_2021_HUMAN ")
            
        return path
    
    def run(self, option, config):
        if( self._validate_input(config) ):
            if( option in range(4) ):
                with open(config, 'r') as g:
                    cfg = json.load(g)
                
                for e in cfg['experiments']:
                    ide = e['identifier']
                    print(f'Experiment [{ide}]')
                    
                    normalization, gtf = self._check_norm(e)
                    geneset = self._check_geneset(e)
                    
                    flagop = ( option==0 ) 
                    if( flagop or option==1):
                        print('\tStep 1 - Running data processing')
                        a = ProcessPathwayScores( e['folder'], e['expression_file'], e['identifier'], normalization, gtf, geneset )
                        a.run()
                        
                    if( flagop or option==2): # drug-pathway-gene- ScoringMatrix
                        flag = self._validate_labels(e)
                        if(flag):
                            a = BuildScoringMatrix( e['folder'], e['identifier'], e['labels_file'], geneset )
                            a.run()
                        
                    if( flagop or option==3): # Training and drug ranking
                        flag = self._validate_labels(e)
                        if(flag):
                            a = DrugRankingAnalysis( e['folder'], e['identifier'], e['labels_file'] )
                            a.run()
                
            else:
                print('Error - Invalid option')
            
import sys

#op = int(sys.argv[1])
#config = sys.argv[2]

import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description=' Caliscoma - Pipeline for drug ranking based on computed pathway scores of disease and healthy samples', formatter_class=RawTextHelpFormatter)

parser.add_argument("-cf", "--configuration_file", action="store", help="(For both modes) Folder to store the files (use the folder where the required files can be found, ex.: /home/user/experiment/ )\n")

parser.add_argument("-rt", "--running_step", action="store", help="0 - Run all steps\n\
1 - Run step 1: Data processing\n\
1 - Run step 2: Model training & modified pathway score matrix\n\
2 - Run step 3: Drug ranking", type=int)

args = parser.parse_args()

a = Pipeline_drugResponseCalibration()
a.run( args.running_step, args.configuration_file)
