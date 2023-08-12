import optuna
import pandas as pd

def compute_scoring_matrix( approved_drugs, w1, w2, w3):
    print("\t\tComputing modified pathway scores for disease samples")
    
    lbfile = 'data_icgc/lirijp_sample_labels.tsv'
    fout = 'data_icgc/lirijp/'
    ide='lirijp'
    
    df = pd.read_csv( f'{fout}/{ide}_pathway_scores.tsv', sep='\t')
    lb = pd.read_csv(lbfile, sep='\t')
    disease = lb[ lb['label']==1 ]['sample'].unique()
    dsa = df[ df['Name'].isin(disease) ][ ['Name', 'Term', 'ES'] ]
    
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
    gb = dsa.groupby(['Name', 'Term'])
    samples = dsa['Name'].unique()
    
    drugs = set(rel_drug_path.keys())
    drugs = drugs.intersection( set(approved_drugs) )
    
    f = open( 'temp_scores.tsv', 'w' )
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
    
    df = pd.read_csv( 'temp_scores.tsv', sep='\t' )
    
    return df

def evaluate(mcs, trained_model):
    gb = mcs.groupby( ['drug', 'sample'] )
    drugs = set( df['drug'].unique() )
    samples = set( df['sample'].unique() )
    
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
        ratio = changed / len(samples)
        
        obj = { 'drug': drug, 'label_changed_ratio': ratio }
        resdf = pd.concat([resdf, pd.DataFrame([ obj ])], ignore_index=True)
    
    score = len( resdf[ resdf['label_changed_ratio'] > 0.7 ] ) / len( resdf )
    return score

def objective(trial):
    fout = 'data_icgc/lirijp/'
    
    trained_model = joblib.load( f'{fout}/selected_model')
    
    # approved drugs for liver cancer - https://www.cancer.gov/about-cancer/treatment/drugs/liver
    ds="Atezolizumab,Avastin (Bevacizumab),Bevacizumab,Cabometyx (Cabozantinib-S-Malate),Cabozantinib-S-Malate,Cyramza (Ramucirumab),Durvalumab,Futibatinib,Imfinzi (Durvalumab),Imjudo (Tremelimumab-actl),Infigratinib Phosphate,Ipilimumab,Keytruda (Pembrolizumab),Lenvatinib Mesylate,Lenvima (Lenvatinib Mesylate),Lytgobi (Futibatinib),Nexavar (Sorafenib Tosylate),Nivolumab,Opdivo (Nivolumab),Pemazyre (Pemigatinib),Pembrolizumab,Pemigatinib,Ramucirumab,Regorafenib,Sorafenib Tosylate,Stivarga (Regorafenib),Tecentriq (Atezolizumab),Tremelimumab-actl,Truseltiq (Infigratinib Phosphate),Yervoy (Ipilimumab)"
    dst = set( [ x.lower().split('-')[0].split(' ')[0] if x.find(' (')==-1 else x.split(' (')[1].replace(')','').split('-')[0].split(' ')[0].lower() for x in ds.split(',') ] )
    dgb = pd.read_csv( 'workflow/drugbank_mapping.tsv', sep='\t')
    approved_drugs = dgb[ dgb['name'].str.lower().isin(dst) ]['id'].values
    
    w1 = trial.suggest_int("w1", 0, 30, log=True)
    w2 = trial.suggest_int("w2", 0, 30, log=True)
    w3 = trial.suggest_int("w3", 0, 30, log=True)

    mcs = compute_scoring_matrix( approved_drugs, w1, w2, w3)
    score = evaluate(mcs, trained_model)
    
    return score

def run_optimization():
    # example: https://github.com/optuna/optuna-examples/blob/main/sklearn/sklearn_simple.py
    
    study = optuna.create_study(direction="maximize")
    study.optimize(objective, n_trials=100)
    print(study.best_trial)

#compute_scoring_matrix( 20, 5, 10)
run_optimization()
