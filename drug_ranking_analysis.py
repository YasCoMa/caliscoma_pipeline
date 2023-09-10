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
    

class DrugRankingAnalysis:
    
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
        self.dsa = df[ ['Name', 'Term', 'ES'] ]
        self.samples = set( self.dsa['Name'].unique() )
        self.dspathways = set( self.dsa['Term'].unique() )
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
                     
    def _get_combined_drug_pathway_score(self):
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
        df = pd.DataFrame( columns=[ 'drug', 'pathway', 'mean_score' ] )
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
                df = pd.concat([df, pd.DataFrame([ obj ])], ignore_index=True)
                rel[d][ obj['pathway'] ] = obj['mean_score']
            i+=1
                
        df.to_csv( f'{fout}/drug_pathway_score_table.tsv', sep='\t', index=None )   
        
        return rel
        
    def compute_scoring_matrix(self, w1, w2, w3):
        print("\t\tComputing modified pathway scores for disease samples")
        
        fout = self.folder_out
        
        rel_drug_path = self._get_combined_drug_pathway_score()
        self._get_pathway_transfer_mapping()
        
        dfmean = pd.read_csv( f'{self.tmeans}', sep='\t', index_col=0 )[ ['abs_diff_mean'] ]
        gb = self.dsa.groupby(['Name', 'Term'])
        samples = self.dsa['Name'].unique()
        
        drugs = list(rel_drug_path.keys())
        
        f = open( f'{fout}/modified_disease_score_by_drug.tsv', 'w' )
        f.write('drug\tpathway\tsample\tmes\n')
        for drug in tqdm( drugs ):
            for pathway in rel_drug_path[ drug ]:
                drug_pathway_score = rel_drug_path[ drug ][ pathway ]
                signal = 1
                if( drug_pathway_score < 0 ):
                    signal = -1
                elif( drug_pathway_score == 0 ):
                    signal = 0
                
                diff_mean_pathway = dfmean.loc[ self.mapping_transfer[pathway], 'abs_diff_mean']
                
                for sample in samples:
                    sample_original_score = gb.get_group( (sample, pathway) )['ES'].values[0]
                    
                    modified_score = sample_original_score
                    
                    if( drug_pathway_score != 0 ):
                        modified_score *= signal
                        
                        if( diff_mean_pathway > 0):
                            q2 = np.quantile( dfmean['abs_diff_mean'], 0.5 )
                            q3 = np.quantile( dfmean['abs_diff_mean'], 0.75 )
                            
                            if( diff_mean_pathway > q3 ):
                                modified_score *= w1
                            elif( diff_mean_pathway >= q2 ):
                                modified_score *= w2
                            else:
                                modified_score *= w3
                            
                    f.write( "%s\t%s\t%s\t%.6f\n" %(drug, pathway, sample, modified_score) )
        f.close()
            
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
    
    def run(self, weights):
        if(self.flag):
            self.compute_scoring_matrix( weights[0], weights[1], weights[2] )
            self.evaluate_rank_drug()
