import os
import optuna
import joblib
import numpy as np
import pandas as pd
import gseapy as gp
from tqdm import tqdm

workflow_path = os.environ.get('path_workflow')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

class WeightOptimization:
    
    def __init__(self, folder, identifier, label_file, drug_list_file, geneset):
        self.flag = True
        
        if(folder[-1]=='/'):
            folder = folder[:-1]
        self.folder = folder
            
        self.id = identifier
        
        self.geneset = geneset
        
        self.folder_out = folder+'/'+identifier
        
        self.label_file = label_file
        
        self.drugs = open( f'{self.folder_out}/{drug_list_file}', "r").read().split("\n")
        
        lbfile = f'{folder}/{self.label_file}'
        lb = pd.read_csv(lbfile, sep='\t')

    def _fill_original_score(self, sample):
        grp = self.grouped_original_scores.get_group(sample)
        paths = grp['Term'].values
        scores = grp['ES'].values
        obj={}
        for p,s in zip(paths, scores):
            obj[p] = s
            
        return obj

    def _get_means(self):
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
        self.dspathways = set( self.dsa['Term'].unique() )
        self.grouped_original_scores = self.dsa.groupby('Name')
        
        hmean = self.hsa[ ['Term', 'ES'] ].groupby('Term').mean()
        dmean = self.dsa[ ['Term', 'ES'] ].groupby('Term').mean()
        diff = dmean['ES'] - hmean['ES']
        
        if( not os.path.isfile( f'{fout}/table_means.tsv' ) ):
            dff = pd.DataFrame( index=hmean.index )
            dff['diff_mean'] = diff
            dff['abs_diff_mean'] = [ abs(x) for x in dff['diff_mean'].values ]
            dff['healthy_mean'] = hmean['ES']
            dff['disease_mean'] = dmean['ES']
            dff.to_csv( f'{fout}/table_means.tsv', sep='\t' )
        self.tmeans = f'{fout}/table_means.tsv'
                     
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
        
        if( not os.path.isfile(f'{fout}/drug_pathway_score_table.tsv') ):
            df.to_csv( f'{fout}/drug_pathway_score_table.tsv', sep='\t', index=None )  
        
        return rel
        
    def compute_scoring_matrix(self, approved_drugs, w1, w2, w3):
        fout = self.folder_out
        
        print("\t\tComputing modified pathway scores for disease samples")
        
        rel_drug_path = {}
        if( not os.path.isfile(f'{fout}/drug_pathway_score_table.tsv') ):
            rel_drug_path = self._get_combined_drug_pathway_score()
        else:
            df = pd.read_csv( f'{fout}/drug_pathway_score_table.tsv', sep='\t')
            for index, row in df.iterrows():
                d = row['drug']
                p = row['pathway']
                s = row['mean_score']
                if(not d in rel_drug_path):
                    rel_drug_path[d]={}
                rel_drug_path[d][p]=s
        
        dfmean = pd.read_csv( f'{self.tmeans}', sep='\t', index_col=0 )[ ['abs_diff_mean'] ]
        #gb = dsa.groupby(['Name', 'Term'])
        gb = self.dsa.groupby('Term')
        samples = self.dsa['Name'].unique()
        ns = len(samples)
        
        drugs = set(rel_drug_path.keys())
        drugs = drugs.intersection( set(approved_drugs) )
        
        df = pd.DataFrame( columns = ['drug', 'pathway', 'sample', 'mes'] )
        for drug in tqdm( drugs ):
            for pathway in rel_drug_path[ drug ]:
                drug_pathway_score = rel_drug_path[ drug ][ pathway ]
                signal = 1
                if( drug_pathway_score < 0 ):
                    signal = -1
                elif( drug_pathway_score == 0 ):
                    signal = 0
                    
                diff_mean_pathway = dfmean.loc[pathway, 'abs_diff_mean']
                q2 = np.quantile( dfmean['abs_diff_mean'], 0.5 )
                q3 = np.quantile( dfmean['abs_diff_mean'], 0.75 )
                
                factor = 1
                if( drug_pathway_score != 0 ):
                    factor *= signal
                    
                    if( diff_mean_pathway > q3 ):
                        factor *= w1
                    elif( diff_mean_pathway >= q2 ):
                        factor *= w2
                    else:
                        factor *= w3
                        
                scores = list( map( lambda x: x*factor, gb.get_group( pathway )['ES'].values ) )
                samples = gb.get_group( pathway )['Name'].values
                temp = pd.DataFrame()
                temp['drug'] = [drug]*ns
                temp['pathway'] = [pathway]*ns
                temp['sample'] = samples
                temp['mes'] = scores
                df = pd.concat([ df, temp ])
                
        return df

    def evaluate(self, df, trained_model):
        fout = self.folder_out
        
        gb = df.groupby( ['drug', 'sample'] )
        drugs = set( df['drug'].unique() )
        samples = set( df['sample'].unique() )
        
        resdf = pd.DataFrame(columns=['drug','label_changed_ratio'])
        for drug in tqdm( drugs ):
            xtest = []
            for sample in samples:
                row = self._fill_original_score(sample)
                svals = [ str(x) for x in list(row.values()) ]
                log = f'{sample}\tbefore\t'+('\t'.join( svals ))+'\n'
                
                paths = gb.get_group( (drug, sample) )['pathway'].values
                scores = gb.get_group( (drug, sample) )['mes'].values 
                for p,s in zip(paths, scores):
                    row[ p ] = s
                xtest.append( list( row.values() ) )
                
                svals = [ str(x) for x in list(row.values()) ]
                log += f'{sample}\tafter\t'+('\t'.join( svals ))+'\n'
                with open( f'{fout}/temp_trials_samples_scores.txt', 'a') as f:
                    f.write(log)
            
            predictions = trained_model.predict( xtest )
            changed = 0
            for s,p in zip(samples, predictions):
                if(p < 0.2):
                    changed+=1
            ratio = changed / len(samples)
            
            obj = { 'drug': drug, 'label_changed_ratio': ratio }
            resdf = pd.concat([resdf, pd.DataFrame([ obj ])], ignore_index=True)
        
        score = len( resdf[ resdf['label_changed_ratio'] > 0.7 ] ) / len( resdf )
        
        print('---> score:', score )
        return score

    def objective(self, trial):
        fout = self.folder_out
        
        trained_model = joblib.load( f'{fout}/selected_model')
        
        f=open( f'{fout}/temp_trials_samples_scores.txt','w')
        f.close()
        
        ct_drugs = self.drugs
        
        w1 = trial.suggest_int("w1", 1, 30)
        w2 = trial.suggest_int("w2", 1, 30)
        w3 = trial.suggest_int("w3", 1, 30)

        #mcs = compute_scoring_matrix( approved_drugs, w1, w2, w3)
        mcs = self.compute_scoring_matrix( ct_drugs, w1, w2, w3)
        score = self.evaluate(mcs, trained_model)
        
        return score

    def run(self):
        fout = self.folder_out
        
        self._get_means();
        
        study = optuna.create_study(direction="maximize")
        study.optimize( self.objective, n_trials=100)
        print(study.best_trial)
        
        ws = study.best_trial.params
        f=open( f"{fout}/best_weights.tsv", "w")
        f.write("weight\tvalue\n")
        for k in ws:
            f.write("%s\t%.6f\n" %(k, ws[k]) )
        f.close()
        

#compute_scoring_matrix( 20, 5, 10)
#run_optimization()
