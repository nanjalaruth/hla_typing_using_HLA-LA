nextflow.enable.dsl=2

params.input = "/cbio/projects/013/custom-bam.ruth/selected/rest/*/*.bam"
params.reference_genome = "/users/nanje/miniconda3/opt/hla-la/graphs/PRG_MHC_GRCh38_withIMGT"
graph_ch = Channel.fromPath(params.reference_genome)
params.outdir = "./output"

process hla_typing {
    tag "Performing HLA typing using HLA-LA"
    label "medium"
    cache "lenient"
    
    input:
        tuple val(dataset), path(reads), path(index), path(graph)

    output:
        path("${dataset}_symlink")

    script:
        out = "/users/kir-luo/ypz679/devel/hla-la_working_dir"
        hla_perl_folder = "/users/kir-luo/ypz679/devel/HLA-LA/src"
        // out = "/scratch3/users/nanje/HLA-LA"

        """
        #HLA typing script
        ${hla_perl_folder}/HLA-LA.pl --BAM ${reads} --graph ${graph} --sampleID ${dataset} --workingDir ${out} --maxThreads 10
        ln -s ${out}/${dataset} ${dataset}_symlink
        """
        
}

process createHpedFiles {
    tag "creating hped files"
    cache "lenient"

    input:
	tuple val(dataset), path(best_guess_folder)

    output:
	path("${best_guess_folder}.ped")

    script:
        """
        #Extract column 3 
        cut -f 3 ${best_guess_folder}/hla/R1_bestguess_G.txt > ${dataset}.test
        #Transpose column to row
        tr "\\n" "\\t" < ${dataset}.test > ${dataset}.test1
        #remove 1st column
        cut -f 2-17 ${dataset}.test1 > ${dataset}.ped
        #add 0 for the first 4 columns
        for i in {1..4}; do sed -i "s/^/0\\t/" ${dataset}.ped ; done
        #add the sample name column twice
        for i in {1..2}; do sed -i "s/^/${dataset}\\t/" ${dataset}.ped ; done 
        #add population name to the end of the file
        sed -i "s/\$/\\tGGVP/" ${dataset}.ped
        """    

}

process concatenateHpedFiles{
    publishDir "${params.outdir}/hla_types", mode: 'copy', overwrite: true
    tag "concatenating hped files"
    cache "lenient"

    input:
	path ped_files

    output:
	path "GGVP.hped"

    script:
	"""
	cat *.ped | \
    sed '1 i FID\\tIID\\tPID\\tMID\\tSEX\\tPHENO\\tA.1\\tA.2\\tB.1\\tB.2\\tC.1\\tC.2\\tDQA1.1\\tDQA1.2\\tDQB1.1\\tDQB1.2\\tDRB1.1\\tDRB1.2\\tDPA1.1\\tDPA1.2\\tDPB1.1\\tDPB1.2\\tPop' > GGVP.hped
	"""
}    

process createCoverageTable {
    tag "creating coverage table"
    cache "lenient"

    input:
	tuple val(dataset), path(best_guess_folder)

    output:
	path("${best_guess_folder}.coverage")

    script:
        """
        #Extract column 3 
        cut -f 3,6 ${best_guess_folder}/hla/R1_bestguess_G.txt > ${dataset}.test
        #Transpose column to row
        tr "\\n" "\\t" < ${dataset}.test > ${dataset}.test1
        #remove 1st column
        cut -f 3- ${dataset}.test1 > ${dataset}.coverage
        #add the sample name column twice
        sed -i "s/^/${dataset}\\t/" ${dataset}.coverage
        """    

}

process concatenateCoverageFiles{
    publishDir "${params.outdir}/hla_types", mode: 'copy', overwrite: true
    tag "concatenating coverage files"
    cache "lenient"

    input:
	path coverage_files

    output:
	path "GGVP.coverage"

    script:
	"""
	cat *.coverage | \
    sed '1 i IID\\tA.1\\tcoverage\\tA.2\\tcoverage\\tB.1\\tcoverage\\tB.2\\tcoverage\\tC.1\\tcoverage\\tC.2\\tcoverage\\tDQA1.1\\tcoverage\\tDQA1.2\\tcoverage\\tDQB1.1\\tcoverage\\tDQB1.2\\tcoverage\\tDRB1.1\\tcoverage\\tDRB1.2\\tcoverage\\tDPA1.1\\tcoverage\\tDPA1.2\\tcoverage\\tDPB1.1\\tcoverage\\tDPB1.2\\tcoverage' > GGVP.coverage
	"""
}  

workflow{    
   //BAM files 
    // input_ch = Channel.fromPath([params.input])
    //        .map{bam -> [file(bam).getSimpleName(), file(bam), file("${bam}.bai")]}
    //        .combine(graph_ch)
    // input_ch.view()

    // Type HLA alleles
    // CRAM files
    input_ch = Channel.fromPath([params.input])
        .map{bam -> [file(bam).getSimpleName(), file(bam), file("${bam}.crai")]}
        .combine(graph_ch)
    out_ch = hla_typing(input_ch)
    
    // Create HLA ped files
    ped_ch = out_ch
        .map{folder -> [file(folder).getSimpleName(), file(folder)]}
    conc_ch = createHpedFiles(ped_ch)

    // Concatenate ped files
    hped_files = conc_ch.collect()
    concatenateHpedFiles(hped_files)

    // Create Coverage files
    cov_ch = out_ch
        .map{folder -> [file(folder).getSimpleName(), file(folder)]}
    conca_ch = createCoverageTable(cov_ch)

    // Concatenate coverage files
    cov_files = conca_ch.collect()
    concatenateCoverageFiles(cov_files)
}
