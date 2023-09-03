import os
import json
import joblib
import numpy as np
import pandas as pd
import gseapy as gp
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score
from tqdm import tqdm

workflow_path = os.environ.get('path_workflow')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

class BuildScoringMatrix:
    
    def __init__(self, folder, identifier, label_file, geneset, model, tmeans):
        self.flag = True
        
        if(folder[-1]=='/'):
            folder = folder[:-1]
        self.folder = folder
            
        self.id = identifier
        
        self.folder_out = folder+'/'+identifier
        
        self.label_file = label_file
        
        self.model = model
        self.tmeans = tmeans
        
        self.geneset = geneset
        
        if( not os.path.isfile( f'{self.folder_out}/{identifier}_pathway_scores.tsv' ) ):
            self.flag = False
            print (f'Error - {identifier}: Pathway score file was not found. You have to run the previous step of the pipeline.')
            
    def _prepare_dataset_xy(self):
        X = []
        y = []
        
        gb = self.hsa.groupby(['Name'])
        for sample in gb.indices.keys():
            X.append( gb.get_group(sample)['ES'].values )
            y.append(0)
        
        gb = self.dsa.groupby(['Name'])
        for sample in gb.indices.keys():
            X.append( gb.get_group(sample)['ES'].values )
            y.append(1)
            
        return X, y

    def _get_means_split_samples(self):
        folder = self.folder
        fout = self.folder_out
        ide = self.id
        lbfile = f'{folder}/{self.label_file}'
        
        df = pd.read_csv( f'{fout}/{ide}_pathway_scores.tsv', sep='\t')
        lb = pd.read_csv(lbfile, sep='\t')
        
        healthy = lb[ lb['label']==0 ]['sample'].unique()
        disease = lb[ lb['label']==1 ]['sample'].unique()
        self.hsa = df[ df['Name'].isin(healthy) ][ ['Name', 'Term', 'ES'] ]
        self.dsa = df[ df['Name'].isin(disease) ][ ['Name', 'Term', 'ES'] ]
        
        X, y = self._prepare_dataset_xy( )
        
        if( not os.path.isfile( f'{fout}/table_means.tsv' ) ):
            hmean = self.hsa[ ['Term', 'ES'] ].groupby('Term').mean()
            dmean = self.dsa[ ['Term', 'ES'] ].groupby('Term').mean()
            diff = dmean['ES'] - hmean['ES']
            
            dff = pd.DataFrame( index=hmean.index )
            dff['diff_mean'] = diff
            dff['abs_diff_mean'] = [ abs(x) for x in dff['diff_mean'].values ]
            dff['healthy_mean'] = hmean['ES']
            dff['disease_mean'] = dmean['ES']
            dff.to_csv( f'{fout}/table_means.tsv', sep='\t' )
            self.tmeans = f'{fout}/table_means.tsv'
            
        return X, y
        
    def train_model(self):
        print("\t\tTraining model")
        
        fout = self.folder_out
        
        X, y = self._get_means_split_samples()

        skf = StratifiedKFold(n_splits=10, shuffle=True)
        iterator = tqdm( skf.split(X, y), desc='Training model to classify disease from healthy samples')
        
        metrics = { 'accuracy': [], 'roc_auc': [], 'precision': [], 'recall': [] }
        
        models = []
        ma=0
        mai=0
        ind = 0
        for train_indexes, test_indexes in iterator:
            X_train = [ X[train_index] for train_index in train_indexes]
            X_test = [ X[test_index] for test_index in test_indexes]
            y_train = [ y[train_index] for train_index in train_indexes]
            y_test = [ y[test_index] for test_index in test_indexes]
            
            m = ElasticNetCV(cv=10, random_state=0)
            m.fit(X_train, y_train)
            y_pred = m.predict(X_test)
            y_pred = [ 1 if x > 0.8 else 0 for x in y_pred ]
            
            for met in metrics:
                metrics[met].append( eval( f"{met}_score(y_test, y_pred)" ) )
            
            models.append( m )
            if( ma < metrics['precision'][-1] ):
                ma = metrics['precision'][-1]
                mai = ind
            ind+=1
            
        trdf = pd.DataFrame()
        for met in metrics:
            trdf[met] = metrics[met]
        trdf.to_csv( f'{fout}/training_evaluation_results.tsv', sep='\t', index=None )
        
        joblib.dump( models[mai], f'{fout}/selected_model')
        self.model =  f'{fout}/selected_model'
                     
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
            
        dspathways = set( self.hsa['Term'].values ).union( self.dsa['Term'].values )
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
            for p in dspathways:
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
                    
                diff_mean_pathway = dfmean.loc[pathway, 'abs_diff_mean']
                
                for sample in samples:
                    sample_original_score = gb.get_group( (sample, pathway) )['ES'].values[0]
                    
                    modified_score = sample_original_score
                    
                    if( drug_pathway_score != 0 ):
                        modified_score *= signal
                        
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
    
    def run(self, weights):
        if(self.flag):
            if( self.label_file!=None ):
                self.train_model()
            self.compute_scoring_matrix( weights[0], weights[1], weights[2] )
            
                    
        
