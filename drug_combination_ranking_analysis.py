import os
import json
import joblib
import numpy as np
import pandas as pd
import gseapy as gp
from tqdm import tqdm

workflow_path = os.environ.get('path_workflow')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]
    
class DrugCombinationAnalysis:
    
    def __init__(self, folder, identifier, label_file, model, tmeans, n_features_model, geneset):
        self.flag = True
        
        self.geneset = geneset
        
        if(folder[-1]=='/'):
            folder = folder[:-1]
        self.folder = folder
            
        self.id = identifier
        
        self.folder_out = folder+'/'+identifier
        
        self.label_file = label_file
        
        if( model == None ):
            model = f"{self.folder_out}/selected_model"
        self.model = model
        
        if( tmeans == None ):
            tmeans = f"{self.folder_out}/table_means.tsv"
        self.tmeans = tmeans
        self.n_features_model = n_features_model
        
        df = pd.read_csv( f'{self.folder_out}/{identifier}_pathway_scores.tsv', sep='\t')
        if( label_file!=None ):
            lbfile = f'{folder}/{self.label_file}'
            lb = pd.read_csv(lbfile, sep='\t')
            disease = lb[ lb['label']==1 ]['sample'].unique()
            df = df[ df['Name'].isin(disease) ]
        df = df[ ['Name', 'Term', 'ES'] ]
        self.samples = set( df['Name'].unique() )
        self.dspathways = set( df['Term'].unique() )
        self.grouped_original_scores = df.groupby('Name')
        dsa=None
        df=None
    
    def _get_pathway_transfer_mapping(self):
        gmt_lib = self.geneset
        gmt = gp.parser.download_library(gmt_lib, 'Human')
        
        mp = {}
        dfmean = pd.read_csv( f'{self.tmeans}', sep='\t', index_col=0 )
        original_paths = set( dfmean.index )
        for ptarget in self.dspathways:
            gt = 0
            chosen = ptarget
            if( not ptarget in original_paths ):
                for pmodel in original_paths:
                    interlen = len( set( gmt[ptarget] ).intersection(set( gmt[pmodel] )) )
                    if( interlen > gt ):
                        gt = interlen
                        chosen = pmodel
            mp[ptarget] = chosen
        self.mapping_transfer = mp
                     
    def _load_drug_pathway_score(self):
        fout = self.folder_out
        
        df = pd.read_csv( f'{workflow_path}/filtered_relation_gene_drug.tsv', sep='\t' )
        mp = {}
        for index, row in df.iterrows():
            drug = row['drug']
            target = row['target']
            if(not drug in mp):
                mp[ drug ] = {}
            mp[ drug ][ target ] = row['relation']
            
        gmt_lib = self.geneset
        gmt = gp.parser.download_library(gmt_lib, 'Human')
        
        rel = {}
        total = len(mp.keys())
        i=1
        for d in mp:
            #print(i, '/', total)
            rel[d] = {}
            valid = []
            for p in self.dspathways:
                targets = mp[d].keys()
                found_pathway_targets = list( set( targets ).intersection( gmt[p] ) )
                
                if( len( found_pathway_targets ) > 0 ) :
                    list_scores = list( map( lambda x: mp[d][x], found_pathway_targets ))
                    mean_score = sum(list_scores) / len(found_pathway_targets)
                    obj = { 'drug': d, 'pathway': p, 'mean_score': mean_score }
                    valid.append(obj)
                    
            for obj in valid:
                rel[d][ obj['pathway'] ] = obj['mean_score']
            i+=1
                
        self.drug_relation_score = rel
        
    def _load_table_means(self):
        self.dfmean = pd.read_csv( f'{self.tmeans}', sep='\t', index_col=0 )[ ['abs_diff_mean'] ]
    
    def _get_drugScore_factor(self, drug, pathway):
        signal = 1
        if( drug in self.drug_relation_score ):
            if( pathway in self.drug_relation_score[ drug ] ):
                drug_pathway_score = self.drug_relation_score[ drug ][ pathway ]
                signal = 1
                if( drug_pathway_score < 0 ):
                    signal = -1
                elif( drug_pathway_score == 0 ):
                    signal = 0
        return signal
    
    def _get_tmeans_factor(self, pathway):
        factor = 1
        pathway = self.mapping_transfer[pathway]
        if( pathway in self.dfmean.index ):
            diff_mean_pathway = self.dfmean.loc[pathway, 'abs_diff_mean']
            q2 = np.quantile( self.dfmean['abs_diff_mean'], 0.5 )
            q3 = np.quantile( self.dfmean['abs_diff_mean'], 0.75 )
            
            if( diff_mean_pathway > q3 ):
                factor *= self.weights[0]
            elif( diff_mean_pathway >= q2 ):
                factor *= self.weights[1]
            else:
                factor *= self.weights[2]
        return factor
        
    def _fill_original_score(self, sample):
        grp = self.grouped_original_scores.get_group(sample)
        paths = grp['Term'].values
        scores = grp['ES'].values
        obj={}
        for p,s in zip(paths, scores):
            obj[p] = s
            
        return obj
    
    def modify_score(self, sample, drug_list):
        row = self._fill_original_score(sample)
        
        for pathway in row:
            for drug in drug_list:
                signal = self._get_drugScore_factor(drug, pathway)
                weight = self._get_tmeans_factor(pathway)
                row[pathway] *= signal*weight
        return row
    
    def _load_drug_combinations(self, drug_comb_list):
        drugs = open( drug_comb_list, 'r' ).read().split('\n')
        drugs = list( filter( lambda x: x!='', drugs ) )
        drugs = list( map(lambda x: x.replace(' ',''), drugs ) )
        return drugs
        
    def evaluate_rank_drug_combination(self, drug_comb_list, w1, w2, w3):
        self.weights = [w1, w2, w3]
        self._load_table_means()
        self._load_drug_pathway_score()
        self._get_pathway_transfer_mapping()
        
        ide = self.id
        fout = self.folder_out
        samples = self.samples
        trained_model = joblib.load( f'{self.model}')
        drugs = self._load_drug_combinations(drug_comb_list)
        
        print("\t\tEvaluating whether the drug-based modified pathways changed the sample condition state")
        
        prlist = {}
        for sample in samples:
            prlist[sample] = {}
        
        resdf = pd.DataFrame(columns=['drug','label_changed_ratio'])
        for drug in tqdm( drugs ):
            xtest = []
            for sample in samples:
                items = drug.split(',')
                row = self.modify_score(sample, items)
                    
                temp = list( row.values() ) 
                if( self.label_file == None ):
                    if( len(temp) > self.n_features_model):
                        temp = temp[:self.n_features_model]
                    if( len(temp) < self.n_features_model):
                        while len(temp) < self.n_features_model:
                            temp.append(0)
                        
                xtest.append( temp )
            
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
        resdf.to_csv( f'{fout}/ranked_drug_combination_evaluation.tsv', sep='\t', index=None )
        
        # Prioritization specific for samples
        df = pd.read_csv( f'{workflow_path}/drugbank_mapping.tsv', sep='\t')
        mp = {}
        for index, row in df.iterrows():
            mp[ row['id'] ] = row['name']
        
        drugs = resdf[ resdf['label_changed_ratio'] > 0.7 ]['drug'].values
        f=open(  f'{fout}/drug_combination_prioritization_list_samples.tsv','w')
        ndr = [ d if( not d in mp ) else f"{d} ({mp[d]})" for d in drugs ]
        f.write( 'sample\t'+( '\t'.join( ndr ) )+'\n' )
        for s in samples:
            vals = [ '1' if prlist[s][x] < 0.2 else '0' for x in drugs ]
            f.write( s+'\t'+( '\t'.join( vals ) )+'\n' )
        f.close()
    
    def run(self, drug_comb_list, weights):
        if(self.flag):
            self.evaluate_rank_drug_combination( drug_comb_list, weights[0], weights[1], weights[2] )
