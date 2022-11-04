nextflow.enable.dsl=2

params.input = "/cbio/projects/013/custom-bam.ruth/selected/rest/*/*.bam"
// params.graph = "/users/nanje/miniconda3/opt/hla-la/graphs/PRG_MHC_GRCh38_withIMGT"

// graph_ch = Channel.fromPath(params.graph)

//params.input = "/cbio/projects/013/custom-bam.ruth/selected/*/*.bam"
// params.index = "/cbio/projects/013/custom-bam.ruth/selected/*/*.bam.bai"
params.reference_genome = "/users/nanje/miniconda3/opt/hla-la/graphs/PRG_MHC_GRCh38_withIMGT"
graph_ch = Channel.fromPath(params.reference_genome)


process hla_typing {
    tag "Performing HLA typing using HLA-LA"
    publishDir "./output/hla_types", mode: 'copy', overwrite: true
    label "medium"
    
    input:
        tuple val(dataset), path(reads), path(index), path(graph)

    script:
        out = "/users/kir-luo/ypz679/devel/hla-la_working_dir"
        hla_perl_folder = "/users/kir-luo/ypz679/devel/HLA-LA/src"

        """
        #HLA typing script
        ${hla_perl_folder}/HLA-LA.pl --BAM ${reads} --graph ${graph} --sampleID ${dataset} --workingDir ${out} --maxThreads 10

        #Extract column 3 
        #cut -f 3 ${out}/${dataset}/hla/R1_bestguess_G.txt > test
        #Transpose column to row
        #tr "\\n" "\\t" < test > test1
        #remove 1st column
        #cut -f 2-17 test1 > test2
        #add 0 for the first 4 columns
        #for i in {1..4}; do sed -i 's/^/0\\t/' test2 ; done
        #add the sample name column twice
        #for i in {1..2}; do sed -i 's/^/${dataset}\\t/' test2 ; done 
        # Combine all the output
        #cat test2 > GGVP.hped
        """        
}

workflow{

    // input_ch = Channel
    //     .fromPath( params.input )
    //     .map { row -> [ file(row).baseName, [ file( row) ] ] }
    // .combine(graph_ch)
    
   //BAM files 
    //input_ch = Channel.fromPath([params.input])
      //      .map{bam -> [file(bam).getSimpleName(), file(bam), file("${bam}.bai")]}
        //    .combine(graph_ch)
    // input_ch.view()

    //CRAM files
    input_ch = Channel.fromPath([params.input])
        .map{bam -> [file(bam).getSimpleName(), file(bam), file("${bam}.crai")]}
        .combine(graph_ch)
    // input_ch.view()

    hla_typing(input_ch)
}
