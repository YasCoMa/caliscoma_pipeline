import os
import json
import numpy as np
import pandas as pd

from rnanorm import UQ
from rnanorm import FPKM
from rnanorm import TPM
from rnanorm import TMM
from rnanorm import CPM
from rnanorm.annotation import GTF

from gseapy import Biomart
import gseapy as gp

workflow_path = os.environ.get('path_workflow')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

class ProcessPathwayScores:
    
    def __init__(self, folder, expression_file, identifier, normalization, gtf, geneset):
        self.flag = True
        
        if(folder[-1]=='/'):
            folder = folder[:-1]
            
        self.id = identifier
        
        self.folder = folder
        
        self.normalization = normalization
        self.gtf = gtf
        self.geneset = geneset
        
        self.folder_out = folder+'/'+identifier
        if( expression_file.endswith('.tsv.gz') ):
            self.raw_data = pd.read_csv( f'{folder}/{expression_file}', compression='gzip', sep='\t'  )
        elif( expression_file.endswith('.tsv') ):
            self.raw_data = pd.read_csv( f'{folder}/{expression_file}', sep='\t'  )
        else:
            self.flag=False
            print('Error - The expression file could not be loaded, check if it is a .tsv or tsv.gz (compressed) file.')
    
        if(self.flag):
            cols = self.raw_data.columns.tolist()
            if( not ( ('submitted_sample_id' in cols) and ('gene_id' in cols) and ('raw_read_count' in cols) ) ):                
                self.flag=False
                print('Error - The mandatory columns were not found in the expression file')
            else:
                if( not os.path.isdir( self.folder_out ) ):
                    os.system( f"mkdir {self.folder_out}" ) 
        
    def _get_gene_names(self):
        df = self.raw_data
        folder = self.folder
        fout = self.folder_out
        ide = self.id

        if(not os.path.isfile( f'{fout}/{ide}_genes.json') ):
            genes = [ x.split('.')[0] for x in list(df['gene_id'].unique()) ]
            with open( f'{fout}/{ide}_genes.json','w') as f:
                json.dump(genes, f)
        else:
            with open( f'{fout}/{ide}_genes.json','r') as f:
                genes = json.load(f)
                
        bm = Biomart()
        
        mapp={}
        if(not os.path.isfile( f'{workflow_path}/mapp_ids.tsv') ):
            f=open( f'{workflow_path}/mapp_ids.tsv','w')
            f.write('ensembl\thgnc_symbol\n')
            f.close()
        else:
            f=open(f'{workflow_path}/mapp_ids.tsv','r')
            for line in f:
                l=line.split('\t')
                if( l[0] != 'ensembl' ):
                    mapp[ l[0] ] = l[1]
            f.close()
            
        genes = list( set(genes) - set(mapp.keys()) )
        
        offset = 200
        i=0
        while (i < len(genes) ):
            print('Stage ', i, '/', len(genes) )
            gs = genes[ i:(i+offset) ]
            queries = {'ensembl_gene_id': gs } # need to be a dict object
            results = bm.query(dataset='hsapiens_gene_ensembl', attributes=['ensembl_gene_id', 'external_gene_name' ], filters=queries)
            for j in results.index:
                id_ = results.iloc[j, 0]
                val = results.iloc[j, 1]
                if(not id_ in mapp):
                    mapp[ id_ ] = val
                    with open( f'{workflow_path}/mapp_ids.tsv','a') as f:
                        f.write("%s\t%s\n" %(id_, val) )
            i+=offset
                
        #print( 'Mapped: ', len(mapp.keys()), '/', len(genes) )

    def _test_map_genes(self, df):
        genes = df['gene_id'].values
        
        ens = 0
        for g in genes:
            if( g.startswith('ENSG') ):
                ens+=1
                
        if( ens>0 ):
            df['gene_id'] = df['gene_id'].apply( lambda x: x.split('.')[0] ) 
            self._get_gene_names()
            
            m = pd.read_csv( f'{workflow_path}/mapp_ids.tsv', sep='\t')
            m = m[ ~m['hgnc_symbol'].isna() ]
            fm = m.groupby('hgnc_symbol').count()
            fm = fm[ fm['ensembl']==1 ]
            m = m[ m['hgnc_symbol'].isin(fm.index) ]
            passed = m['ensembl'].values
            m = m.reset_index()
            mp = {}
            for index, row in m.iterrows():
                mp[ row['ensembl'] ] = row['hgnc_symbol']

            df = df[ df['gene_id'].isin(passed) ]
            df['gene_id'] = [ mp[x] for x in list(df['gene_id']) ]
        return df

    def get_gct_rnanorm_files(self):
        print("\t\tGenerating input for rnanorm and the gct file")
    
        df = self.raw_data
        folder = self.folder
        fout = self.folder_out
        ide = self.id
        
        df = df[ ['gene_id','submitted_sample_id','raw_read_count' ] ]
        df = self._test_map_genes(df)
        
        self.raw_data = None
        
        fil = df.groupby('gene_id').count().sort_values(by='raw_read_count')
        fil = fil[ fil['submitted_sample_id'] <= len(df['submitted_sample_id'].unique()) ]
        n = fil.tail(1).iloc[0,0]
        indexok = list( fil[ fil['raw_read_count'] > (n/2) ].index ) # removing genes that are not present in at least half of the sample ids
        df = pd.pivot_table( df, values='raw_read_count', index=['gene_id'], columns=['submitted_sample_id'], aggfunc = np.sum) # generating gct table for single sample analysis in gseapy
        df = df.loc[indexok, :].fillna(0)

        df.reset_index(inplace=True)
        cols = df.columns.tolist()
        cols[0]='Name'
        df.columns = cols
        df['Description'] = df['Name']
        df = df[ ['Name','Description']+cols[1:] ]
        df = df[ df['Name'].isin(indexok) ].fillna(0)

        # input rnorm lirijp
        df = df.transpose()
        head = df.values[0,:]
        index = list(df.index)[2:]
        df = pd.DataFrame(df.values[2:,:], columns=head)
        df.index = index
        df.to_csv( f'{fout}/input_rnanorm.csv')
        
    def _check_gene_length(self, genes):
        gtf_file = f'{workflow_path}/gtf_{self.gtf}.gtf.gz'
        gtflen = GTF(gtf=gtf_file).length
        
        mapp = pd.read_csv( f'{workflow_path}/mapp_ids.tsv', sep='\t')
        mp={}
        for i in mapp.index:
            mp[ mapp.iloc[i, 0] ] = mapp.iloc[i, 1]
        filt = list( filter( lambda x: x.split('.')[0] in mp.keys(),  gtflen.keys() ))
        
        newgl = {}
        for g in filt:
            hgnc = mp[ g.split('.')[0] ]
            if(hgnc in genes):
                newgl[ hgnc ] = gtflen[g]
        se = pd.Series(newgl)
        genesok = list( newgl.keys() )
        
        return se, genesok
        
    def normalize_counts(self):
        print("\t\tGenerating normalized table for ssGSEA")
        
        folder = self.folder
        fout = self.folder_out
        ide = self.id
        
        df = pd.read_csv( f'{fout}/input_rnanorm.csv', index_col=0)
        df = df.fillna(0)
        
        if( self.normalization in ['fpkm_uq', 'cpm', 'tmm'] ):
            if(self.normalization == 'fpkm_uq'):
                transformer = UQ()
            if(self.normalization == 'cpm'):
                transformer = CPM()
            if(self.normalization == 'tmm'):
                transformer = TMM()
        else:
            genes = df.columns
            gene_lengths, genesok = self._check_gene_length( genes )
            df = df[ genesok ]
            
            if(self.normalization == 'tpm'):
                transformer = TPM( gene_lengths = gene_lengths )
            elif(self.normalization == 'fpkm'):
                transformer = FPKM( gene_lengths = gene_lengths )
                
        transformer.fit(df)
        index = df.index
        cols = df.columns
        df = transformer.transform(df)
        df = pd.DataFrame(df, columns = cols, index = index)
        
        df = df.transpose()
        df.reset_index(inplace=True)
        cols = df.columns.tolist()
        cols[0]='NAME'
        df.columns=cols
        df['Description'] = df['NAME']
        df = df[ ['NAME', 'Description']+cols[1:] ]
        df.to_csv( f'{fout}/table_normalized.gct', sep='\t', index=None)

    def run_ssgsea(self):
        print("\t\tRunning ssGSEA analysis to obtain the pathway scores")
        
        folder = self.folder
        fout = self.folder_out
        ide = self.id
        
        gmt_lib = self.geneset
        
        gmt = gp.parser.download_library(gmt_lib, 'Human')
        genesok = set()
        for k in gmt.keys():
            genesok = genesok.union( set(gmt[k]) )
        
        df = pd.read_csv( f'{fout}/table_normalized.gct', sep='\t')
        cols = df.columns.tolist()[2:]
        df['Gene'] = [ x.upper() for x in df['Description'].values ]
        df['NAME'] = [ x.upper() for x in df['Description'].values ]
        df = df[ ['Gene', 'NAME']+cols ]
        
        genes = set(df['Gene'])
        passed = genes.intersection(genesok)
        df = df[ df['NAME'].isin(passed) ]
        df.to_csv( f'{fout}/final_{ide}_table.tsv', sep='\t', index=None)
        df = None
        
        ss = gp.ssgsea(data = f'{fout}/final_{ide}_table.tsv', gene_sets=gmt, outdir=None, sample_norm_method='rank', no_plot=True)
        df = ss.res2d
        df.to_csv( f'{fout}/{ide}_pathway_scores.tsv', sep='\t', index=None)
    
    def run(self):
        if(self.flag):
            self.get_gct_rnanorm_files()
            self.normalize_counts()
            self.run_ssgsea()
        
