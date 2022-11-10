nextflow.enable.dsl=2

params.input = "/cbio/projects/013/custom-bam.ruth/selected/rest/*/*.bam"
params.reference_genome = "/users/nanje/miniconda3/opt/hla-la/graphs/PRG_MHC_GRCh38_withIMGT"
graph_ch = Channel.fromPath(params.reference_genome)

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
        """    

}

process concatenateHpedFiles{
    publishDir "./output/hla_types", mode: 'copy', overwrite: true
    tag "concatenating hped files"
    cache "lenient"

    input:
	path ped_files

    output:
	path "GGVP.hped"

    script:
	"""
	cat *.ped | \
    sed '1 i FID\tIID\tPID\tMID\tSEX\tPHENO\tA.1\tA.2\tB.1\tB.2\tC.1\tC.2\tDQA1.1\tDQA1.2\tDQB1.1\tDQB1.2\tDRB1.1\tDRB1.2\tDPA1.1\tDPA1.2\tDPB1.1\tDPB1.2' > GGVP.hped
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
    
    // Concatenate ped files
    conc_ch = createHpedFiles(ped_ch)
    hped_files = conc_ch.collect()
    concatenateHpedFiles(hped_files)
}
