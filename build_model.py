import os
import json
import joblib
import numpy as np
import pandas as pd

from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score
from tqdm import tqdm

workflow_path = os.environ.get('path_workflow')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

class BuildModel:
    
    def __init__(self, folder, identifier, label_file):
        self.flag = True
        
        if(folder[-1]=='/'):
            folder = folder[:-1]
        self.folder = folder
            
        self.id = identifier
        
        self.folder_out = folder+'/'+identifier
        
        self.label_file = label_file
        
        if( not os.path.isfile( f'{self.folder_out}/{identifier}_pathway_scores.tsv' ) ):
            self.flag = False
            print (f'Error - {identifier}: Pathway score file was not found. You have to run the previous step of the pipeline.')
        else:
            df = pd.read_csv( f'{self.folder_out}/{identifier}_pathway_scores.tsv', sep='\t')
            self.dsa = df[ ['Name', 'Term', 'ES'] ]
            self.samples = set( self.dsa['Name'].unique() )
            self.dspathways = set( self.dsa['Term'].unique() )
            
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
            diff = hmean['ES'] - dmean['ES']
            
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
    
    def run(self):
        if(self.flag):
            self.train_model()
        
