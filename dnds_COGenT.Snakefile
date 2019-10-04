"""
This snakemake pipeline takes as input a singleton hmmsearch result and outputs a
phylogeny (Raxml) built from an alignment of ribosomal proteins (mafft).
This technique was used to create a ribosomal reference tree 
in Hehemann et al., 2016 (doi:10.1038/ncomms12860).

"""

from Bio import SeqIO
import glob
import pandas as pd


configfile:
    "dnds_COGenT_config.yml"

contig_wildcard_string = config["protein_directory"] + "/*" + config["unaligned_nt_ext"]
prots = [f.split('/')[-1].replace(config["unaligned_nt_ext"], "") for f in glob.glob(contig_wildcard_string)]
rule target:    
    input:
        config["output_directory"] + 'dnds_consolidated.csv'
        #config["output_directory"] + config['result_file_prefix'] + ".dnds.csv"

rule consolidate_results:
    input:
        expand(config["output_directory"] + "{prot}.dnds.csv", prot=prots)
    output: 
        config["output_directory"] + 'dnds_consolidated.csv'
    run:
        all_results = [pd.read_csv(f) for f in input]
        df = pd.concat(all_results)
        df.to_csv(str(output), index=False)


rule get_dn_ds:
    input:
        config["output_directory"] + "{prot}.pal2nal.nt.fasta"
        #config["output_directory"] + config['result_file_prefix'] + ".pal2nal.nt.fasta"
    output:
        config["output_directory"] + "{prot}.dnds.csv"
    shell:
        "Rscript dnds_COGenT.calc.R {input} {config[group_file]} {config[group_type]} {output}"

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