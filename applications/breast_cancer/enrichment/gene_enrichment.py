# Gene set over-representation analysis
"""
Gene sets downloaded from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H

- c2.all.v2025.1.Hs.symbols.gmt (C2 curated gene sets)
	- Description: Gene sets in this collection are curated from various sources, including 
		online pathway databases and the biomedical literature. Many sets are also contributed 
		by individual domain experts. The gene set page for each gene set lists its source. The 
		C2 collection is divided into the following two subcollections: Chemical and genetic 
		perturbations (CGP) and Canonical pathways (CP)
"""

import gseapy as gp
import pandas as pd
import re

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 0)
pd.set_option('display.expand_frame_repr', False)


data_path = '../data/cleaned_data/cleaned_data.csv'
df = pd.read_csv(data_path)

target_names = ['histological_type', 'pathologic_stage', 'status']
feature_names = df.columns.tolist()

#----------------------------------------------------------------
# Genes
#----------------------------------------------------------------
gene_pathways = {
	# histologic subtype
	'cdh1': ['CDH1', 'CES2', 'ABHD1', 'CES8', 'KIAA0174', 'AP1G1', 'USP10',
		'DYNC1LI2', 'CTCF', 'ATXN1L', 'ATMIN', 'AIP', 'TERF2IP', 'GABARAPL2'],
	'snrpn': ['SNRPN', 'IPW', 'SNORD116-20', 'SNORD116-28', 'PAR-SN', 'SNURF', 'PAR5'],

	# stage
	'hook3': ['HOOK3', 'QKI', 'SGK196', 'MYST3', 'MAP4K5', 'FNTA', 'HGSNAT','AGPAT6', 'GOLGA7', 'IKBKB', 
		'VDAC3', 'WHSC1L1', 'SLC20A2'],

	# status
	'vps72': ['VPS72', 'SF3B4', 'MRPL9', 'MRPL24', 'PSMB4', 'PSMD4'],
	'loc162632': ['LOC162632', 'FAM106A', 'FAM106C', 'CCDC144A', 'CCDC144B', 'LOC220594'],
	'ankle2': ['ANKLE2', 'GOLGA3', 'EP400', 'POLE', 'SETD1B', 'MLL2', 'C12orf51', 'ULK1', 'SFRS8'],

	# miRNA modules (gene targets)
	'mir133a1': ['FBXO2', 'ACTC1', 'DES', 'MYH11', 'LOC728264'],
	'mir210': ['P4HA1', 'ISCU', 'RIC8A', 'CA9', 'EPS8L2', 'BTNL9', 'NDRG1', 'CD151', 'PDDC1', 'PKP3', 
		'PHRF1', 'BET1L', 'NAP1L4']
}

#----------------------------------------------------------------
# Over-representation analysis
#----------------------------------------------------------------

pathway_name = 'hook3'

# background set (all genes in analyzed dataset)
background_set = sorted({f.replace('*','') for f in feature_names if f.endswith('*') and f not in target_names })

selected_genes = sorted(set(gene_pathways[pathway_name]))
print(f'Pathway: {pathway_name}')
print(f'Size of background set: {len(background_set)}')
print(f'Pathway genes: {selected_genes}')

# Use only stable gene set servers
gene_sets = '../data/gene_sets/'

gene_sets += 'c2.all.v2025.1.Hs.symbols.gmt'
# gene_sets += 'c2.cp.v2025.1.Hs.symbols.gmt'
# gene_sets += 'c2.cp.reactome.v2025.1.Hs.symbols.gmt'
# gene_sets += 'c3.all.v2025.1.Hs.symbols.gmt'
# gene_sets += 'c3.mir.v2025.1.Hs.symbols.gmt'
# gene_sets += 'c5.all.v2025.1.Hs.symbols.gmt'
# gene_sets += 'c5.go.bp.v2025.1.Hs.symbols.gmt'
# gene_sets += 'c6.all.v2025.1.Hs.symbols.gmt'

print(f'Gene set: {gene_sets}')
print(f'----------------------------------------------------------------')

enr = gp.enrich(gene_list=selected_genes, gene_sets=gene_sets, background=background_set, outdir=None)

# check if no results returned
if not isinstance(enr.results, pd.DataFrame) or enr.results.empty:
	print('No enrichment results returned.')
	exit()

res = enr.results[['Term', 'Overlap', 'Adjusted P-value', 'Genes']]

# only consider gene sets where the minimum overlap with selected genes is min_overlap
min_overlap = 3
mask = res['Overlap'].str.split('/').str[0].astype(int) >= min_overlap
res = res.loc[mask]
res = res.sort_values('Adjusted P-value')

print(res.head(5).to_string(index=False))

