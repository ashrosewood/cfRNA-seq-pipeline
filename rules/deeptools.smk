rule index_TE:
    input:
        "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam"
    output:
        "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam.bai"
    conda:
        "../envs/omic_qc_wf.yaml"
    shell:
        """samtools index -@ 4 {input} {output}"""

#rule cpm_tracks:
#    input:
#        bam = "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam",
#        idx = "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam.bai"
#    output:
#        wig = "samples/bigwig_TE/{sample}_cpm.bw"
#    params:
#        exclude = config["blackList"]
#    conda:
#        "../envs/deeptools.yaml"
#    shell:
#        "bamCoverage -p 4 --normalizeUsing CPM -bs 1 -b {input.bam} -o {output.wig}"

rule fwd_tracks:
    input:
        bam = "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam",
        idx = "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
       wig = "samples/bigwig_TE/{sample}_fwd.bw"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -p 4 --normalizeUsing CPM --filterRNAstrand forward -bs 1 -b {input.bam} -o {output.wig}"

rule rev_tracks:
    input:
        bam = "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam",
        idx = "samples/star_TE/{sample}/Aligned.sortedByCoord.out.bam.bai"
    output:
        wig = "samples/bigwig_TE/{sample}_rev.bw"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -p 4 --normalizeUsing CPM --filterRNAstrand reverse -bs 1 -b {input.bam} -o {output.wig}"

rule uniq_fwd_tracks:
    input:
        bam = "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        idx = "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    output:
       wig = "samples/star/{sample}_bam/{sample}_fwd.bw"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -p 4 --normalizeUsing CPM --filterRNAstrand forward -bs 1 -b {input.bam} -o {output.wig}"

rule uniq_rev_tracks:
    input:
        bam = "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam",
        idx = "samples/star/{sample}_bam/Aligned.sortedByCoord.out.bam.bai"
    output:
        wig = "samples/star/{sample}_bam/{sample}_rev.bw"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -p 4 --normalizeUsing CPM --filterRNAstrand reverse -bs 1 -b {input.bam} -o {output.wig}"

