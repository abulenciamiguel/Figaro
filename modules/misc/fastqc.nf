process fastqcRawSE {
        container 'staphb/fastqc:0.12.1'

        tag "Check quality of raw ${sample}"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq)

        output:
        tuple val(sample), path("*fastqc*"), emit: qualRaw


        script:
        """
        fastqc --outDir . $fastq
        """
}


process fastqcTrimmedSE {
        container 'staphb/fastqc:0.12.1'

        tag "Check quality of trimmed ${sample}"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq)

        output:
        tuple val(sample), path("*fastqc*"), emit: qualTrimmed


        script:
        """
        fastqc --outDir . $fastq
        """
}




process fastqcRawPE {
        container 'staphb/fastqc:0.12.1'

        tag "Checking quality of ${sample}"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq_1), path(fastq_2)

        output:
        tuple val(sample), path("*fastqc*"), emit: qualRaw


        script:
        """
        fastqc --outDir . $fastq_1 $fastq_2
        """
}


process fastqcTrimmedPE {
        container 'staphb/fastqc:0.12.1'

        tag "Checking quality of ${sample}"

        publishDir (
        path: "${params.outDir}/${task.process.replaceAll(":","_")}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(fastq_1), path(fastq_2)

        output:
        tuple val(sample), path("*fastqc*"), emit: qualTrimmed


        script:
        """
        fastqc --outDir . $fastq_1 $fastq_2
        """
}