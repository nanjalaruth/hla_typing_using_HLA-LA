nextflow.enable.dsl=2

// params.input = "/cbio/projects/013/custom-bam.ruth/selected/*/*.bam"
// params.index = "/cbio/projects/013/custom-bam.ruth/selected/*/*.bam.bai"
params.graph = "/users/nanje/miniconda3/opt/hla-la/graphs/PRG_MHC_GRCh38_withIMGT"

// index_ch = Channel.fromPath(params.index)
graph_ch = Channel.fromPath(params.graph)

process hla_typing {
    tag "Performing HLA typing using HLA-LA"
    publishDir "./output", mode: 'copy', overwrite: false
    label "medium"
    
    input:
        tuple val(dataset), path(reads), path(index), path(graph)
    output:
        path(output)
    script:
        out = "/scratch3/users/nanje/HLA-LA/output"
        """
        HLA-LA.pl --BAM ${reads} --graph ${graph} --sampleID ${dataset} --workingDir ${out} --maxThreads 10
        """        
}

workflow{
    // input_ch = Channel
    //     .fromPath( params.input )
    //     .map { row -> [ file(row).baseName, [ file( row) ] ] }
    // .combine(graph_ch)
   input_ch = Channel.fromPath(["/cbio/projects/013/custom-bam.ruth/selected/*/*.bam"])
        .map{bam -> [file(bam).getSimpleName(), file(bam), file("${bam}.bai")]}
        .combine(graph_ch)

    hla_typing(input_ch)
}
