import optuna
import joblib
import numpy as np
import pandas as pd
from tqdm import tqdm

class WeightOptimization:
    
    def __init__(self, folder, identifier, label_file, drug_list_file):
        self.flag = True
        
        if(folder[-1]=='/'):
            folder = folder[:-1]
        self.folder = folder
            
        self.id = identifier
        
        self.folder_out = folder+'/'+identifier
        
        self.label_file = label_file
        
        self.drugs = open( f'{self.folder_out}/{drug_list_file}', "r").read().split("\n")
        
        lbfile = f'{folder}/{self.label_file}'
        lb = pd.read_csv(lbfile, sep='\t')
        
        disease = lb[ lb['label']==1 ]['sample'].unique()
        dsa = df[ df['Name'].isin(disease) ][ ['Name', 'Term', 'ES'] ]
        self.grouped_original_scores = df.groupby('Name')

    def _fill_original_score(sample):
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
        
    def compute_scoring_matrix( approved_drugs, w1, w2, w3):
        fout = self.folder_out
        
        print("\t\tComputing modified pathway scores for disease samples")
        
        df = pd.read_csv( f'{fout}/drug_pathway_score_table.tsv', sep='\t')
        rel_drug_path = {}
        for index, row in df.iterrows():
            d = row['drug']
            p = row['pathway']
            s = row['mean_score']
            if(not d in rel_drug_path):
                rel_drug_path[d]={}
            rel_drug_path[d][p]=s
        
        dfmean = pd.read_csv( f'{fout}/table_means.tsv', sep='\t', index_col=0 )[ ['abs_diff_mean'] ]
        #gb = dsa.groupby(['Name', 'Term'])
        gb = dsa.groupby('Term')
        samples = dsa['Name'].unique()
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

    def evaluate(df, trained_model):
        gb = df.groupby( ['drug', 'sample'] )
        drugs = set( df['drug'].unique() )
        samples = set( df['sample'].unique() )
        
        resdf = pd.DataFrame(columns=['drug','label_changed_ratio'])
        for drug in tqdm( drugs ):
            xtest = []
            for sample in samples:
                row = _fill_original_score(sample)
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
        
        print('---> score:', score, resdf['label_changed_ratio'].values )
        return score

    def objective(trial):
        fout = self.folder_out
        
        trained_model = joblib.load( f'{fout}/selected_model')
        
        f=open( f'{fout}/temp_trials_samples_scores.txt','w')
        f.close()
        
        ct_drugs = self.drugs
        
        w1 = trial.suggest_int("w1", 1, 30)
        w2 = trial.suggest_int("w2", 1, 30)
        w3 = trial.suggest_int("w3", 1, 30)

        #mcs = compute_scoring_matrix( approved_drugs, w1, w2, w3)
        mcs = compute_scoring_matrix( ct_drugs, w1, w2, w3)
        score = evaluate(mcs, trained_model)
        
        return score

    def run_optimization():
        fout = self.folder_out
        
        if( not os.path.isfile( f'{fout}/table_means.tsv' ) ):
            self._get_means();
        
        study = optuna.create_study(direction="maximize")
        study.optimize(objective, n_trials=100)
        print(study.best_trial)
        
        ws = study.best_trial
        f=open( f"{fout}/best_weights.tsv", "w")
        f.write("weight\tvalue\n")
        for k in ws:
            f.write("%s\t%.6f\n" %(k, ws[k]) )
        f.close()
        

#compute_scoring_matrix( 20, 5, 10)
#run_optimization()
