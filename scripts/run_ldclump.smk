# load config
configfile: "scripts/config.yaml"

# Get the gene list
with open(config["gene_input"], "r") as f:
    genes = [line.strip() for line in f if line.strip()]

# For each gene, run the clumping
rule run_all:
	input:
		expand("{outdir}_thresh{ldthresh}/clumped_all.txt", outdir=config["outdir"], ldthresh=config["ldthresh"])
    
rule run_LD_clump:
    output:
        "{outdir}_thresh{ldthresh}/per_gene/{gene}_clumped_format.txt" 
    params:
        plink_prefix=config["plink_prefix"],
        plink_suffix=config["plink_suffix"],
        outdir=config["outdir"],
        ldthresh=config["ldthresh"],
        gene_variant_input=config["gene_variant_input"]
    threads: 1
    singularity:
        "/software/hgi/softpack/installs/groups/macromapsqtl//macromapsqtl/4-scripts/singularity.sif"
    shell:
        r"""
            # Make out dir
            mkdir -p {params.outdir}_thresh{params.ldthresh}
            mkdir -p {params.outdir}_thresh{params.ldthresh}/per_gene
            mkdir -p temp

            # Extract the variants of interest from the gene x variant list. Put in temporary file
            awk -v g="{wildcards.gene}" '$1 == g' {params.gene_variant_input} | awk '{{print $2}}' | uniq > temp/{wildcards.gene}_variants.txt

            # Check how many variants, if >1 clump
            num_lines=$(wc -l < temp/{wildcards.gene}_variants.txt)

            if [ "$num_lines" -gt 1 ]; then

                echo ".. ~~~~~~~ " $num_lines "variants for this gene, LD clumping ~~~~~~~"
                echo "..Getting chromosome number"
                num_chr=$(awk -F':' 'NR==1 {{gsub("chr", "", $1); print $1}}' temp/{wildcards.gene}_variants.txt)

                echo "..Filtering plink files"
                plink --bfile {params.plink_prefix}${{num_chr}}{params.plink_suffix} --extract temp/{wildcards.gene}_variants.txt --make-bed --out temp/{wildcards.gene}

                # Get the specific sumstats for this gene
                echo "..Getting specific gene x variant pairs"
                awk -v g="{wildcards.gene}" 'NR==1 || $1 == g' {params.gene_variant_input} >> temp/{wildcards.gene}_assoc.txt

                # LD clump
                echo "..Clumping"
                plink --bfile temp/{wildcards.gene} --clump temp/{wildcards.gene}_assoc.txt --clump-p1 1 --clump-p2 1 --clump-r2 {params.ldthresh} --clump-kb 1000 --out temp/{wildcards.gene}

                # Simplify the output
                echo "..Reformatting"
                Rscript scripts/simplify_clump.r temp/{wildcards.gene}.clumped {params.outdir}_thresh{params.ldthresh}/per_gene/{wildcards.gene}_clumped_format.txt

            else

                echo ".. ~~~~~~~ 1 variants for this gene, NOT LD clumping ~~~~~~~"
                variant=$(awk '{{print $1}}' temp/{wildcards.gene}_variants.txt | head -n 1)
                echo -e "{wildcards.gene}\t${{variant}}\t${{variant}}" > {params.outdir}_thresh{params.ldthresh}/per_gene/{wildcards.gene}_clumped_format.txt

            fi

            # Remove the the intermediate files
            rm temp/{wildcards.gene}*

        """

def gather_per_gene_outputs(wildcards):
    return expand("{outdir}_thresh{ldthresh}/per_gene/{gene}_clumped_format.txt",
                  outdir=config["outdir"], ldthresh=config["ldthresh"], gene=genes)

rule aggregate:
    input:
        gather_per_gene_outputs
    output:
        "{outdir}_thresh{ldthresh}/clumped_all.txt"
    params:
        outdir=config["outdir"],
        ldthresh=config["ldthresh"]
    threads: 1
    shell:
        r"""
            # Add header
            echo -e "phenotype_id\tvariant_id\tqtl_clump_index" > {output}

            # Loop through input
            for f in {params.outdir}_thresh{params.ldthresh}/per_gene/*; do
                cat $f >> {output}
            done
        """