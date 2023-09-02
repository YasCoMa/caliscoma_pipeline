import os
import json
import joblib
import numpy as np
import pandas as pd
from tqdm import tqdm

workflow_path = os.environ.get('path_workflow')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]
    

class DrugRankingAnalysis:
    
    def __init__(self, folder, identifier, label_file, model, tmeans):
        self.flag = True
        
        if(folder[-1]=='/'):
            folder = folder[:-1]
        self.folder = folder
            
        self.id = identifier
        
        self.folder_out = folder+'/'+identifier
        
        self.label_file = label_file
        
        self.model = model
        slef.tmeans = tmeans
        
        if( not os.path.isfile( f'{self.folder_out}/modified_disease_score_by_drug.tsv' ) ):
            self.flag = False
            print (f'Error - {identifier}: The modified pathway score file for disease samples was not found. You have to run the previous step of the pipeline.')
        else:
            df = pd.read_csv( f'{self.folder_out}/{identifier}_pathway_scores.tsv', sep='\t')
            if( lbfile!=None ):
                lbfile = f'{folder}/{self.label_file}'
                lb = pd.read_csv(lbfile, sep='\t')
                disease = lb[ lb['label']==1 ]['sample'].unique()
                df = df[ df['Name'].isin(disease) ]
            df = df[ ['Name', 'Term', 'ES'] ]
            self.grouped_original_scores = df.groupby('Name')
            dsa=None
            df=None
            
    def _fill_original_score(self, sample):
        grp = self.grouped_original_scores.get_group(sample)
        paths = grp['Term'].values
        scores = grp['ES'].values
        obj={}
        for p,s in zip(paths, scores):
            obj[p] = s
            
        return obj
    
    def evaluate_rank_drug(self):
        print("\t\tEvaluating whether the drug-based modified pathways changed the sample condition state")
        
        ide = self.id
        fout = self.folder_out
        
        df = pd.read_csv( f'{fout}/modified_disease_score_by_drug.tsv', sep='\t' )
        gb = df.groupby( ['drug', 'sample'] )
        drugs = set( df['drug'].unique() )
        samples = set( df['sample'].unique() )
            
        trained_model = joblib.load( f'{self.model}')
        
        prlist = {}
        for sample in samples:
            prlist[sample] = {}
        
        resdf = pd.DataFrame(columns=['drug','label_changed_ratio'])
        for drug in tqdm( drugs ):
            xtest = []
            for sample in samples:
                row = self._fill_original_score(sample)
                
                paths = gb.get_group( (drug, sample) )['pathway'].values
                scores = gb.get_group( (drug, sample) )['mes'].values 
                for p,s in zip(paths, scores):
                    row[ p ] = s
                xtest.append( list( row.values() ) )
            
            predictions = trained_model.predict( xtest )
            changed = 0
            for s,p in zip(samples, predictions):
                if(p < 0.2):
                    changed+=1
                prlist[s][drug] = p
            ratio = changed / len(samples)
            
            obj = { 'drug': drug, 'label_changed_ratio': ratio }
            resdf = pd.concat([resdf, pd.DataFrame([ obj ])], ignore_index=True)
        
        resdf = resdf.sort_values(by='label_changed_ratio', ascending=False)    
        resdf.to_csv( f'{fout}/ranked_drug_evaluation.tsv', sep='\t', index=None )
        
        # Prioritization specific for samples
        df = pd.read_csv( f'{workflow_path}/drugbank_mapping.tsv', sep='\t')
        mp = {}
        for index, row in df.iterrows():
            mp[ row['id'] ] = row['name']
        
        drugs = resdf[ resdf['label_changed_ratio'] > 0.7 ]['drug'].values
        f=open(  f'{fout}/prioritization_list_samples.tsv','w')
        ndr = [ d if( not d in mp ) else f"{d} ({mp[d]})" for d in drugs ]
        f.write( 'sample\t'+( '\t'.join( ndr ) )+'\n' )
        for s in samples:
            vals = [ '1' if prlist[s][x] < 0.2 else '0' for x in drugs ]
            f.write( s+'\t'+( '\t'.join( vals ) )+'\n' )
        f.close()
    
    def run(self):
        if(self.flag):
            self.evaluate_rank_drug()
