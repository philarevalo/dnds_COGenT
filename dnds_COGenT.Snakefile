"""
This snakemake pipeline takes as input a singleton hmmsearch result and outputs a
phylogeny (Raxml) built from an alignment of ribosomal proteins (mafft).
This technique was used to create a ribosomal reference tree 
in Hehemann et al., 2016 (doi:10.1038/ncomms12860).

"""

from Bio import SeqIO
import glob
import pandas as pd
from itertools import combinations
import os


configfile:
    "dnds_COGenT_config.yml"

contig_wildcard_string = config["protein_directory"] + "/*" + config["unaligned_nt_ext"]
prots = [f.split('/')[-1].replace(config["unaligned_nt_ext"], "") for f in glob.glob(contig_wildcard_string)]
rule target:    
    input:
        config["output_directory"] + 'dnds_consolidated.csv',
        div=config["output_directory"] + 'div_consolidated.csv'
        #config["output_directory"] + config['result_file_prefix'] + ".dnds.csv"

rule consolidate_results:
    input:
        dnds=expand(config["output_directory"] + "{prot}.dnds.csv", prot=prots),
        div=expand(config["output_directory"] + "{prot}.div.csv", prot=prots)
    output: 
        dnds=config["output_directory"] + 'dnds_consolidated.csv',
        div=config["output_directory"] + 'div_consolidated.csv'
    run:
        all_results = [pd.read_csv(f) for f in input['dnds']]
        df = pd.concat(all_results)
        df.to_csv(str(output['dnds']), index=False)

        all_results = [pd.read_csv(f) for f in input['div']]
        df = pd.concat(all_results)
        df.to_csv(str(output['div']), index=False)


rule get_dn_ds:
    input:
        config["output_directory"] + "{prot}.pal2nal.nt.fasta"
        #config["output_directory"] + config['result_file_prefix'] + ".pal2nal.nt.fasta"
    output:
        config["output_directory"] + "{prot}.dnds.csv"
    shell:
        "Rscript dnds_COGenT.calc.R {input} {config[group_file]} {config[group_type]} {output}"

rule get_diversity:
    input:
        aa=config["output_directory"] + '{prot}.aa.aln.fasta',
        nt=config["output_directory"] + "{prot}.pal2nal.nt.fasta"
    output:
        config["output_directory"] + "{prot}.div.csv"
    run:
        populations = pd.read_csv(config['group_file'], index_col='Strain', sep='\t')
        seqs_aa = {s.id: str(s.seq) for s in SeqIO.parse(str(input['aa']), 'fasta')}
        seqs_nt = {s.id: str(s.seq) for s in SeqIO.parse(str(input['nt']), 'fasta')}
        new_rows = []
        prot = os.path.basename(str(input['aa'])).replace('.aa.aln.fasta', '')
        for i1, i2 in combinations(list(seqs_aa.keys()), 2):
            p1 = str(populations.loc[i1, config['group_type']])
            p2 = str(populations.loc[i2, config['group_type']])
            s1 = seqs_aa[i1]
            s2 = seqs_aa[i2]


            diffs_aa = 0
            length_aa = 0
            for b1, b2 in zip(s1, s2):
                if b1 != '-' and b2 != '-':
                    length_aa += 1
                    if b1 != b2:
                        diffs_aa += 1

            s1 = seqs_nt[i1]
            s2 = seqs_nt[i2]
            diffs_nt = 0
            length_nt = 0
            for b1, b2 in zip(s1, s2):
                if b1 != '-' and b2 != '-':
                    length_nt += 1
                    if b1 != b2:
                        diffs_nt += 1


            new_rows.append([i1, i2, int(p1!=p2), '-'.join(sorted([p1, p2])), diffs_aa, length_aa, diffs_aa/length_aa, diffs_nt, length_nt, diffs_nt/length_nt, prot])

        df = pd.DataFrame(new_rows, columns=['Strain1', 'Strain2', 'diff_group', 'pop_pair', 'AA Differences', 'AA Length', 'AA Divergence', 'NT Differences', 'NT Length', 'NT Divergence', 'prot'])
        df.to_csv(str(output), index=False)


rule align_prots:
    input:
        config["protein_directory"] + "/{prot}" + config["unaligned_aa_ext"]
    output:
        config["output_directory"] + '{prot}.aa.aln.fasta'
    run:
        p = str(output).split('/')[-1].split('_')[0]
        shell("mafft-linsi {input} > {output}")

rule align_prots_nt:
    input:
        nt_filter = config["output_directory"] + "/{prot}" + ".filtered.nt.fasta",
        aa_aln = config["output_directory"] + '{prot}.aa.aln.fasta'
    output:
        config["output_directory"] + "{prot}.pal2nal.nt.fasta"
    run:
        p = str(output).split('/')[-1].split('_')[0]
        shell("perl ./pal2nal.v14/pal2nal.pl {input[aa_aln]} {input[nt_filter]} -output fasta >> {output}")

rule remove_stop_nt:
    input:
        config["protein_directory"] + "/{prot}" + config["unaligned_nt_ext"]
    output:
        config["output_directory"] + "/{prot}" + ".filtered.nt.fasta"
    run:
        new_s = []
        stop_codons = ['TAG', 'TAA', 'TGA']
        for s in SeqIO.parse(str(input), 'fasta'):
            if str(s.seq)[-3:] in stop_codons:
                s.seq = s.seq[0:-3]
            else:
                temp = s
            new_s.append(s)
        SeqIO.write(new_s, str(output), 'fasta')