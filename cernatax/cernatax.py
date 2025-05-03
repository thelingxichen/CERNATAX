import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text


class CERNATAX():
    
    def __init__(self):
                
        wdr = os.path.dirname(os.path.abspath(__file__))        
        self.ref_db = pd.read_csv(os.path.join(wdr, 'db', 'ceRNA_database.csv'), index_col=0)
        
    
    def summarize_ref_db(self):
        ref_db = self.ref_db
        
        print('A total of {} miRNA-mRNA interaction and {} miRNA-lncRNA interaction'.format(
            ref_db[ref_db.type == 'miRNA-mRNA'].shape[0],
            ref_db[ref_db.type == 'miRNA-lncRNA'].shape[0],
        ))
        print(ref_db.type.value_counts())
        
        
    def find_ceRNA_axis_by_DEG(self, deg_df):
        ref_db = self.ref_db
        
        deg_miRNA = deg_df[deg_df.type == 'miRNA'].gene.unique()
        deg_ceRNA = deg_df[deg_df.type != 'miRNA'].gene.unique()
    
        df = ref_db[ref_db.miRNA.isin(deg_miRNA) & ref_db.ceRNA.isin(deg_ceRNA)]
        df = pd.merge(df, deg_df['log2FC'], left_on='miRNA', right_index=True)
        df = pd.merge(df, deg_df['log2FC'], left_on='ceRNA', right_index=True)
        df.columns = ['miRNA', 'ceRNA', 'species', 'database', 'type', 'miRNA_log2FC', 'ceRNA_log2FC']
        df['inference'] = 'strict'
        ceRNA_df = df[(df.miRNA_log2FC * df.ceRNA_log2FC < 0)]
    
        tmp = ceRNA_df.groupby('miRNA').agg({'type': lambda x: len(set(x))})
        axis_df = ceRNA_df[ceRNA_df.miRNA.isin(tmp[tmp.type > 1].index)]
    
        self.ceRNA_df = ceRNA_df
        self.axis_df = axis_df

        return ceRNA_df, axis_df        
    
    def expand_ceRNA_axis_by_loose_DEG(self, deg_strict_df, deg_df):
        # just expand ceRNA based on miRNA axis
        ref_db = self.ref_db
    
        axis_miRNA_list = self.axis_df.miRNA.unique()   
        deg_ceRNA = deg_df[deg_df.type != 'miRNA'].gene.unique()
        
        df = ref_db[ref_db.miRNA.isin(axis_miRNA_list) & ref_db.ceRNA.isin(deg_ceRNA)]
        df = pd.merge(df, deg_strict_df['log2FC'], left_on='miRNA', right_index=True)
        df = pd.merge(df, deg_df['log2FC'], left_on='ceRNA', right_index=True)
        df.columns = ['miRNA', 'ceRNA', 'species', 'database', 'type', 'miRNA_log2FC', 'ceRNA_log2FC']
        df['inference'] = 'loose'
        return df
    
    def plot_cohort_ceRNA_exp(self, adata, groupby, ceRNA_list):
        sc.pl.violin(adata, ceRNA_list, groupby=groupby)

    def cohort_ceRNA_corr(self, adata, ceRNA_axis_list):

        data = []
        for gene1, gene2 in ceRNA_axis_list:
        
            res = stats.pearsonr(adata[:, gene1].X.T[0], adata[:, gene2].X.T[0])
            data.append(['All', gene1, gene2, 'Pearson Correlation', res.statistic, res.pvalue])
    
            res = stats.spearmanr(adata[:, gene1].X.T[0], adata[:, gene2].X.T[0])
            data.append(['All', gene1, gene2, 'Spearman Correlation', res.correlation, res.pvalue])
    
            res = stats.kendalltau(adata[:, gene1].X.T[0], adata[:, gene2].X.T[0])
            data.append(['All', gene1, gene2, "Kendall's tau", res.correlation, res.pvalue])
    
        corr_df = pd.DataFrame(data)
        corr_df.columns = ['patient', 'gene1', 'gene2', 'correlation type', 'correlation', 'p value']
        return corr_df
    
    def screen_GWAS_SNP(self, sig_snp_data, gene_df):
        gene_count = {}
        snp_list = []

        for key, snp_row in sig_snp_data.iterrows():
            SNP_chr = snp_row['CHR']
            SNP_pos = int(snp_row['BP'])
            
            snp_row['gene'] = []
            for gene, row in gene_df.iterrows():
                if SNP_chr != row['chr']:
                    continue
                if SNP_pos < row['start'] or SNP_pos > row['end']:
                    continue
                #print(i, SNP_chr, SNP_pos, gene, row)
                if gene in gene_count:
                    gene_count[gene] += 1
                else:
                    gene_count[gene] = 1
                snp_row['gene'].append(gene)
            if len(snp_row['gene']) > 0:
                snp_list.append(snp_row)
                    
        GWAS_gene_snp_df = pd.DataFrame(snp_list)
        sub_gene_df = gene_df[gene_df.gene.isin(gene_count.keys())]
        sub_gene_df['SNP_count'] = sub_gene_df.gene.apply(lambda x: gene_count[x])
        return GWAS_gene_snp_df, sub_gene_df, gene_count

    def plot_gene_GWAS_SNP(self, GWAS_gene_snp_df, gene_count, ceRNA_list, out_fn):

        fig, axes = plt.subplots(1, 2, figsize=(10, 5))

        for i, gene in enumerate(ceRNA_list):
            ax = axes[i]
            df = GWAS_gene_snp_df[GWAS_gene_snp_df['Gene']==gene]
            sns.scatterplot(df, 
                            x='OR', y='-log10(p)', palette='Set2', #size=1,
                        hue='SNP Type', ax=ax, 
                        hue_order=['A>C', 'A>G', 'C>A', 'C>G', 'C>T', 'G>A', 'T>C', 
                                'Small Insertion', 'Small Deletion'])
            # Calculate y = -log10(0.05)
            threshold = -np.log10(0.05)

            # Add Horizontal Line
            ax.set_title(f'{gene}(n={gene_count[gene]})')
            ax.axhline(y=threshold, color='gray', linestyle='--', linewidth=2, label=f"p = 0.05")
            #ax.xticks(rotation=10)

            # Create text labels
            texts = []
            data = df[df['-log10(p)'] > 2]
            for i in range(len(data)):
                texts.append(ax.text(data['OR'].iloc[i], data['-log10(p)'].iloc[i], data['SNP'].iloc[i], fontsize=10))

            # print(texts)
            # Adjust the positions of the texts to prevent overlap
            adjust_text(texts, ax=ax)


        # Get handles and labels from both plots
        handles1, labels1 = axes[0].get_legend_handles_labels()
        handles2, labels2 = axes[1].get_legend_handles_labels()
        axes[0].legend().remove()
        axes[1].legend().remove()


        # Combine handles and labels, ensuring no duplicates
        handles = handles1 + handles2
        labels = labels1 + labels2
        unique_handles_labels = dict(zip(labels, handles))

        # Add common legend
        fig.legend(unique_handles_labels.values(), unique_handles_labels.keys(), 
                loc='center left', bbox_to_anchor=(0.9, 0.5), ncol=1, frameon=False)
        # Adjust spacing
        #plt.tight_layout()
        fig.savefig(out_fn, format='pdf')
        plt.show()        