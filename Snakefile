import socket
import sys

configfile: "config.yaml"

BPLINK = ["bed", "bim", "fam"]

wildcard_constraints:
  cohort="[^+/]+"

adgc_stem = ('/sc/arion/projects/LOAD/Data/'
             'ADGC/ADGC_2018/Brian_pipeline_withfamily_and_adc8')


rule all:
    input:
        expand('output/PRS_{apoerem}_{prsset}.prsice',
               apoerem=['noAPOE', 'withAPOE'],
               prsset=['superagers', 'nonorthwell']),
        expand('analysis/{sensitivity}/stats/{statfile}',
            statfile = ['fig1stats.txt', 'fig2stats.txt', 'OR_table.tsv'],
            sensitivity = ['agesensitivity', 'primary'])

'''
Filter ADGC phenotypes to new cohorts since Lambert et al. 2013
Output phenotypes and list of samples to keep.
'''

rule filter_ADGC:
  input: adgc_stem + '/ADGC_withAllPCs.withADSPexclusions.pheno.tsv'
  output:
    pheno = 'input/ADGC_phenotypes.tsv',
    ikeep_withIGAP = 'input/withIGAP.ikeep',
    ikeep_withIGAP_qc = 'input/qcpass_withIGAP.ikeep',
    ikeep_noIGAP_qc = 'input/qcpass_noIGAP.ikeep'
  log: 'input/filter_ADGC.log'
  conda: 'rpyenv.yaml'
  script: 'scripts/filter_ADGC.R'


'''
merge and process ADGC and northwell phenotypes
'''

rule prepare_phenos:
    input:
        lz = '../phenotype4brian.txt',
        adgc = rules.filter_ADGC.output.pheno,
        PCA = '/dev/null',
        lz_cases = 'input/cases.irem'
    output:
        phenos = 'merged/phenotypes.tsv',
        hist = 'plots/agehist.pdf'
    log: 'logs/phenotype_prep.log'
    conda: 'rpyenv.yaml'
    script: 'scripts/prepare_phenos.R'


'''
Filter ADGC genotypes to remove QC fails, relateds and IGAP samples
'''

adgc_genos = adgc_stem + '/post_impute/data/all_chrall_filtered_withpheno'

rule filter_ADGC_genos:
    input:
        genos = multiext(adgc_genos, '.bed', '.bim', '.fam'),
        ikeep = rules.filter_ADGC.output.ikeep_noIGAP_qc
    output: multiext('input/ADGC', '.bed', '.bim', '.fam')
    params:
        inplink = adgc_genos,
        outplink = 'input/ADGC'
    conda: 'rpyenv.yaml'
    shell:
        """
plink --keep-allele-order --allow-no-sex \
 --bfile {params.inplink} --keep {input.ikeep} \
 --hwe 1e-10 'midp' --geno 0.05 \
 --make-bed --out {params.outplink}
"""


'''
Make list of NIA-LOAD cases from LZ to exclude.
This is because of potential overlap with ADGC
'''

rule lz_cases:
    input: 'data/demographic_201906.txt'
    output: 'input/cases.irem'
    shell:
        """
awk 'BEGIN {{FS=OFS="\t"}} $1 ~ "^[0-9]{{2}}(M|AD)" {{print $1,$2}}' \
 {input} > {output}
"""


'''
Concatenate the lz chromosomes into a single file
'''

rule lz_concat:
    input:
        expand('data/chr{chrom}.gq_cpra_common.{ext}',
               ext=['bed', 'bim', 'fam'], chrom=list(range(1,23)))
    output:
        mergelist = temp('input/merge_list.txt'),
        genos = temp(multiext('input/chrall.lz', '.bed', '.bim', '.fam'))
    params:
        firstchrom = 'data/chr1.gq_cpra_common',
        otherchroms = expand('data/chr{chrom}.gq_cpra_common',
                             chrom=list(range(2,23))),
        outplink = 'input/chrall.lz'
    conda: 'rpyenv.yaml'
    shell:
        r"""
echo {params.otherchroms} | sed 's/ /\n/g' > {output.mergelist}
plink --keep-allele-order --allow-no-sex --threads 8 --memory 30000 \
 --bfile {params.firstchrom} --merge-list {output.mergelist} \
 --make-bed --out {params.outplink}
"""


'''
Exclude NIA-LOAD cases from LZ genotypes.
This is because of potential overlap with ADGC
'''

rule lz_cases_geno:
    input:
        genos = rules.lz_concat.output.genos,
        irem = rules.lz_cases.output
    output: multiext('input/lz', '.bed', '.bim', '.fam')
    params:
        inplink = rules.lz_concat.params.outplink,
        outplink = 'input/lz'
    conda: 'rpyenv.yaml'
    shell:
        """
plink --keep-allele-order --allow-no-sex \
 --bfile {params.inplink} --remove {input.irem} \
 --make-bed --out {params.outplink}
"""


'''
find overlap in variant positions to improve the
quality and speed of merging
'''

rule overlap:
    input:
        adgc = 'input/ADGC.bim',
        lz = 'input/lz.bim'
    output:
        adgc = 'input/ADGC.overlap',
        lz = 'input/lz.overlap'
    conda: 'rpyenv.yaml'
    script: 'scripts/shared_loc.R'


'''
filter to overlapping positions
'''

rule overlap_filter:
    input:
        overlap = 'input/{cohort}.overlap',
        plink = multiext('input/{cohort}', '.bed', '.bim', '.fam')
    output:
        multiext('input/{cohort}.overlapping', '.bed', '.bim', '.fam')
    params:
        ins = 'input/{cohort}',
        out = 'input/{cohort}.overlapping'
    conda: 'rpyenv.yaml'
    shell:
        """
plink --keep-allele-order --allow-no-sex \
 --bfile {params.ins} --extract {input.overlap} \
 --make-bed --out {params.out}
"""

# MERGE GENOTYPES:

rule Flip:
    input:
        bim = 'input/{cohort}.overlapping.bim',
        bed = 'input/{cohort}.overlapping.bed',
        fam = 'input/{cohort}.overlapping.fam',
        fasta = config["ref"]
    params:
        outdir = 'premerge/{cohort}'
    output:
        temp(expand('premerge/{{cohort}}_flipped.{ext}',
                    ext=BPLINK))
    conda: 'rpyenv.yaml'
    shell:
        """
flippyr -p --plinkMem 20000 \
 -o {params.outdir} {input.fasta} {input.bim}"""

# Recode subcohort plink file to vcf
rule Plink2vcf:
    input: rules.Flip.output
    output: temp('premerge/{cohort}+flipped.vcf.gz')
    params:
        indat = 'premerge/{cohort}_flipped',
        out = 'premerge/{cohort}+flipped'
    shell:
        """
./plink2 --bfile {params.indat} --recode vcf bgz \
--real-ref-alleles --out {params.out}"""

# Index vcf
rule subcohort_Indexvcf:
    input: 'premerge/{cohort}+flipped.vcf.gz'
    output: temp('premerge/{cohort}+flipped.vcf.gz.tbi')
    conda: 'rpyenv.yaml'
    shell: 'bcftools index -t -f {input}'

# Merge
rule merge:
    input:
        vcf = expand('premerge/{cohort}+flipped.vcf.gz',
                     cohort=['lz', 'ADGC']),
        tbi = expand('premerge/{cohort}+flipped.vcf.gz.tbi',
                     cohort=['lz', 'ADGC'])
    output: temp('merged/genos.vcf.gz')
    conda: 'rpyenv.yaml'
    shell:
        r"""
bcftools merge -m none {input.vcf} | \
bcftools annotate --threads 2 --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Oz -o {output}
"""

# Make merged plink files
rule plink:
    input: rules.merge.output
    output: temp(expand('merged/genos_vcfid.{ext}', ext=BPLINK))
    params:
        out_ = 'merged/genos_vcfid'
    shell:
        """
./plink2 --vcf {input} --const-fid --make-bed \
--real-ref-alleles --out {params.out_}"""

rule fam_ref:
    input: expand('input/{cohort}.fam', cohort=['lz', 'ADGC'])
    output: temp('merged/allfam.fam')
    shell: 'cat {input} > {output}'

rule fix_ids:
    input:
        plink = rules.plink.output,
        allfam = rules.fam_ref.output,
        phenos = rules.prepare_phenos.output.phenos
    output: expand('merged/genos.{ext}', ext=BPLINK)
    params:
        ins = rules.plink.params.out_,
        out_ = 'merged/genos',
        status = 'all_status',
        sex = 'sex'
    conda: 'rpyenv.yaml'
    script: 'scripts/fix_fam.R'

rule missingness:
    input: rules.fix_ids.output
    output:
        expand('merged/genos_highcall.{ext}', ext=BPLINK)
    params:
        ins = rules.fix_ids.params.out_,
        out_ = 'merged/genos_highcall'
    conda: 'rpyenv.yaml'
    shell:
        r"""
plink --keep-allele-order --bfile {params.ins} \
 --geno 0.05 --hwe 1e-6 'midp' --make-bed \
 --out {params.out_}
"""

rule missingness_diff:
    input: rules.missingness.output
    output:
        missingrep = 'merged/genos_cr95diff.missing',
        failvars = 'merged/genos_cr95diff.exclude'
    params:
        ins = rules.missingness.params.out_,
        out_ = 'merged/genos_cr95diff'
    conda: 'rpyenv.yaml'
    shell:
        r"""
plink --keep-allele-order --bfile {params.ins} \
 --test-missing 'midp' --out {params.out_}
sed -r 's/[[:blank:]]+/ /g;s/^\s|\s$//g' {output.missingrep} |
  awk '$5 < 1e-4 {{print $2}}' > {output.failvars}
"""

rule missingness_filt:
    input:
        geno = rules.missingness.output,
        excl = rules.missingness_diff.output.failvars
    output: multiext('merged/CR95_noDiffMiss', '.bed', '.bim', '.fam')
    params:
        ins = rules.missingness.params.out_,
        out_ = 'merged/CR95_noDiffMiss'
    conda: 'rpyenv.yaml'
    shell:
        """
plink --keep-allele-order --bfile {params.ins} --exclude {input.excl} \
 --make-bed --out {params.out_}
"""

#PCA stuff below:

# ---- Exlude SNPs with a high missing rate and low MAF----
rule pca_snp_qc:
    input: rules.missingness_filt.output
    output: temp(expand("PCA/SnpQc.{ext}", ext=BPLINK)),
    params:
        stem = rules.missingness_filt.params.out_,
        out = "PCA/SnpQc",
        miss = 0.02,
        MAF = 0.02
    conda: 'rpyenv.yaml'
    shell:
        """
plink --keep-allele-order --bfile {params.stem} \
--geno {params.miss} --maf {params.MAF} \
--make-bed --out {params.out}"""

# ---- Prune SNPs, autosome only ----
#  Pruned SNP list is used for IBD, PCA and heterozigosity calculations

rule PruneDupvar_snps:
    input: rules.pca_snp_qc.output
    output:
        expand("PCA/nodup.{ext}", ext=['prune.in', 'prune.out']),
        "PCA/nodup.dupvar.delete"
    params:
        indat = rules.pca_snp_qc.params.out,
        dupvar = "PCA/nodup.dupvar",
        out = "PCA/nodup"
    conda: 'rpyenv.yaml'
    shell:
        """
plink --keep-allele-order --bfile {params.indat} --autosome --indep 50 5 1.5 \
--list-duplicate-vars --out {params.out}
Rscript scripts/DuplicateVars.R {params.dupvar}"""

# Prune sample dataset
rule sample_prune:
    input:
        rules.pca_snp_qc.output,
        prune = "PCA/nodup.prune.in",
        dupvar = "PCA/nodup.dupvar.delete"
    output:
        temp(expand("PCA/genotypes_pruned.{ext}", ext=BPLINK))
    params:
        indat_plink = rules.pca_snp_qc.params.out,
        out = "PCA/genotypes_pruned"
    conda: 'rpyenv.yaml'
    shell:
        """
plink --keep-allele-order --bfile {params.indat_plink} \
--extract {input.prune} --exclude {input.dupvar} \
--make-bed --out {params.out}"""


# Run PCA to for population stratification
rule PopulationStratification:
    input: rules.sample_prune.output
    output:
        expand("PCA/PCA.{ext}", ext=['eigenval', 'eigenvec'])
    params:
        indat = "PCA/genotypes_pruned",
        out = "PCA/PCA"
    conda: 'rpyenv.yaml'
    shell:
        """
plink --keep-allele-order --bfile {params.indat} --pca 10 \
--out {params.out}
"""

rule updatepheno:
    input:
        pheno = 'merged/phenotypes.tsv',
        eigenvec = "PCA/PCA.eigenvec"
    output: "merged/phenotypes_withJointPCs.tsv"
    conda: 'rpyenv.yaml'
    script: 'scripts/addPCs.R'

'''
Run PRS clumping, thresholding and scoring:
'''

gwasdir = '../APOEe3_full_pipeline_rerun_superaged/partitioned'

gwas = {'noAPOE': gwasdir + '/allstats_1mb_noAPOE.tsv.gz',
        'withAPOE': gwasdir + '/allstats.1mb.tsv.gz'}

prscov = {'noAPOE':  '@jointPC[1-10],sex',
          'withAPOE': '@jointPC[1-10],sex,APOE --cov-factor APOE'}

def pcol(wildcards):
    if wildcards.prsset == 'superagers':
        return 'status'
    return 'status_nonorthwell'

bar_levels = ('5e-08,1e-07,1e-06,1e-05,0.0001,0.001,0.01,0.05,'
              '0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1')

rule do_prs:
    input:
        genotypes = rules.missingness_filt.output,
        phenotypes = 'merged/phenotypes_withJointPCs.tsv',
        gwas = lambda wildcards: gwas[wildcards.apoerem]
    params:
        phenocol = pcol,
        inplink = rules.missingness_filt.params.out_,
        out = 'output/PRS_{apoerem}_{prsset}',
        cov = lambda wildcards: prscov[wildcards.apoerem]
    output: multiext('output/PRS_{apoerem}_{prsset}', '.all.score', '.prsice')
    conda: 'rpyenv.yaml'
    shell:
        """
./PRSice.R \
 --prsice ./PRSice \
 --base {input.gwas} \
 --out {params.out} \
 --target {params.inplink} \
 --cov-file {input.phenotypes} \
 --pheno-file {input.phenotypes} \
 --cov-col {params.cov} \
 --pheno-col {params.phenocol} \
 --beta --binary-target T --stat BETA \
 --bp POS --chr CHR --snp SNP --A1 Allele1 --A2 Allele2 --pvalue P \
 --clump-kb 250 --clump-p 1 --clump-r2 0.1 \
 --model add --perm 10000 \
 --quantile 10 --quant-break 2,4,6,8,10 --quant-ref 6 \
 --fastscore --all-score --bar-levels {bar_levels} \
 --thread 8 --seed 3548238268 --print-snp
"""

'''
combine PRS and phenotypes for analysis
'''

rule process_scores:
    input:
        scores = 'output/PRS_{apoerem}_superagers.all.score',
        phenos = 'merged/phenotypes_withJointPCs.tsv'
    output: 'output/scores+phenos_{apoerem}.Rdata'
    conda: 'rpyenv.yaml'
    script: 'scripts/process_scores_phenos.R'

rule run_analyses:
    input:
        scores = 'output/scores+phenos_noAPOE.Rdata',
        thresh_superager = 'output/PRS_noAPOE_superagers.prsice',
        thresh_standard = 'output/PRS_noAPOE_nonorthwell.prsice',
    params:
        outdir = 'analysis/{sensitivity}'
    output:
        plots = expand('analysis/{{sensitivity}}/plots/{plot}.{ext}',
            plot = ['fig1', 'fig2', 'ageplot', 'lineplot', 'ROCs',
                    'OR', 'OR_nocont', 'OR_strat', 'prsdists'],
            ext = ['png', 'pdf']),
        stats = expand('analysis/{{sensitivity}}/stats/{statfile}',
            statfile = ['fig1stats.txt', 'fig2stats.txt', 'OR_table.tsv'])
    conda: 'rpyenv.yaml'
    script: 'scripts/analysis/main.R'
