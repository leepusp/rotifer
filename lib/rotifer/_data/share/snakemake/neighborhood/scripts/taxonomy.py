from rotifer.db import ncbi
import pandas as pd
input_file = snakemake.input["ipg_combined"]
output_file = snakemake.output["out_file"]

ipg = pd.read_pickle(input_file).drop_duplicates('pid')
a = ncbi.assemblies(taxonomy=False)
a.rename({'#assembly_accession': 'assembly_accession'}, axis=1, inplace=True)
a = a.query('assembly_accession in @ipg.assembly.unique().tolist()')
tc =   ncbi.TaxonomyCursor()
taxonomies = tc.fetchall(a.taxid.tolist())
a = a[['assembly_accession','taxid']]
result = a.merge(taxonomies,
                 how='left').merge(ipg,
                                  left_on='assembly_accession',
                                   right_on='assembly',

                                   how='left'
                                   )
result['position'] = result.nucleotide + ":" + result.start.astype(str) + '-' + result.stop.astype(str)
result.to_pickle(output_file)
