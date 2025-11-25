// enable dsl2
nextflow.enable.dsl=2

// import modules
include {concatenate} from '../modules/ont/concatenate.nf'
include {nanoq} from '../modules/ont/nanoq.nf'
include {fastqcRawSE} from '../modules/misc/fastqc.nf'
include {fastqcTrimmedSE} from '../modules/misc/fastqc.nf'

include {nextcladeDB} from '../modules/misc/nextclade.nf'
include {nextcladeSubtype} from '../modules/misc/nextclade.nf'
include {sierra} from '../modules/misc/sierra.nf'
include {pdfReport} from '../modules/misc/pdfReport.nf'

include {minimap2SE} from '../modules/misc/minimap2.nf'
include {sam2bam} from '../modules/ont/samtools.nf'
include {sortIndexMinimap} from '../modules/ont/samtools.nf'
include {trimPrimer} from '../modules/ont/ivar.nf'
include {sortIndexIvar} from '../modules/ont/samtools.nf'
include {medaka} from '../modules/ont/medaka.nf'


include {getReadNumberLength} from '../modules/misc/runReport.nf'
include {getPhredTrimmed} from '../modules/misc/runReport.nf'
include {getMappedPercent} from '../modules/misc/runReport.nf'
include {getBAMinfo} from '../modules/misc/runReport.nf'
include {htmlRunReport} from '../modules/misc/runReport.nf'



workflow ontAmplicon {
    // Set channel for the fastq directories
    barcodeDirs = file("${params.inDir}/barcode*", type: 'dir', maxdepth: 1 )
    fastqDir = file("${params.inDir}/*.fastq*" , type: 'file', maxdepth: 1)

    main:
        if (barcodeDirs) {
            ch_sample = Channel
                    .fromPath(barcodeDirs)
                    .filter(~/.*barcode[0-9]{1,4}$/)
                    .map{dir ->
                        def reads = 0
                        for (x in dir.listFiles()) {
                            if (x.isFile() && x.toString().contains('.fastq')) {
                                reads += x.countFastq()
                            }
                        }
                        return[dir.baseName, dir]
                    }           
            concatenate(ch_sample)
            fastqcRawSE(concatenate.out.concatFastq)
            nanoq(concatenate.out.concatFastq)
            fastqcTrimmedSE(nanoq.out.trimmedFastq)
            getReadNumberLength(concatenate.out.concatFastq.join(nanoq.out.trimmedFastq))


        } else if (fastqDir) {
            ch_sample = Channel
                    .fromPath(fastqDir)
                    .filter(file -> file.name =~ /.*\.fastq(\.gz)?$/)
                    .map{file ->
                        def baseName = file.name.replaceAll(/\.fastq(\.gz)?$/, '')
                        def reads = file.countFastq()
                        return [baseName, file]
                    }
            fastqcRawSE(ch_sample)
            nanoq(ch_sample)
            fastqcTrimmedSE(nanoq.out.trimmedFastq)
            getReadNumberLength(ch_sample.join(nanoq.out.trimmedFastq))

        } else {
            log.error "Please specify a valid folder containing ONT basecalled, barcoded fastq files or the concatenated fastq files e.g. --inDir ./raw/fastq_pass/ or --inDir ./fastqConcatenated/"
            System.exit(1)
        }

        // Insert here de novo assembly
        // Kraken to get the most likely subtype; Create database containing only HIV sequences
        // Get reference then do reference-based assembly using minimap2
        minimap2SE(nanoq.out.trimmedFastq, params.reference)
        sam2bam(minimap2SE.out.sam)
        sortIndexMinimap(sam2bam.out.bam)
        trimPrimer(sortIndexMinimap.out.bamBai)
        sortIndexIvar(trimPrimer.out.trimmedBam)
        medaka(sortIndexIvar.out.bamBai)

        getPhredTrimmed(nanoq.out.trimmedFastq.join(getReadNumberLength.out.numlenJSON))
        getMappedPercent(minimap2SE.out.sam.join(getPhredTrimmed.out.phredJSON))
        getBAMinfo(sortIndexIvar.out.bamBai.join(getMappedPercent.out.mappedJSON))

        nextcladeDB()
        nextcladeSubtype(medaka.out.consensus, nextcladeDB.out.database)
        sierra(medaka.out.consensus.join(nextcladeSubtype.out.tsv))

        htmlRunReport(getBAMinfo.out.depthJSON.collect(), sierra.out.json.collect{it[0]})
        pdfReport(sierra.out.json, params.markdownFile)
}


workflow ontShotgun {
    // Set channel for the fastq directories
    barcodeDirs = file("${params.inDir}/barcode*", type: 'dir', maxdepth: 1 )
    fastqDir = file("${params.inDir}/*.fastq*" , type: 'file', maxdepth: 1)

    main:
        if (barcodeDirs) {
            ch_sample = Channel
                    .fromPath(barcodeDirs)
                    .filter(~/.*barcode[0-9]{1,4}$/)
                    .map{dir ->
                        def reads = 0
                        for (x in dir.listFiles()) {
                            if (x.isFile() && x.toString().contains('.fastq')) {
                                reads += x.countFastq()
                            }
                        }
                        return[dir.baseName, dir]
                    }           
            concatenate(ch_sample)
            fastqcRawSE(concatenate.out.concatFastq)
            nanoq(concatenate.out.concatFastq)
            fastqcTrimmedSE(nanoq.out.trimmedFastq)
            getReadNumberLength(concatenate.out.concatFastq.join(nanoq.out.trimmedFastq))


        } else if (fastqDir) {
            ch_sample = Channel
                    .fromPath(fastqDir)
                    .filter(file -> file.name =~ /.*\.fastq(\.gz)?$/)
                    .map{file ->
                        def baseName = file.name.replaceAll(/\.fastq(\.gz)?$/, '')
                        def reads = file.countFastq()
                        return [baseName, file]
                    }
            fastqcRawSE(ch_sample)
            nanoq(ch_sample)
            fastqcTrimmedSE(nanoq.out.trimmedFastq)
            getReadNumberLength(ch_sample.join(nanoq.out.trimmedFastq))

        } else {
            log.error "Please specify a valid folder containing ONT basecalled, barcoded fastq files or the concatenated fastq files e.g. --inDir ./raw/fastq_pass/ or --inDir ./fastqConcatenated/"
            System.exit(1)
        }

        minimap2SE(nanoq.out.trimmedFastq, params.reference)
        sam2bam(minimap2SE.out.sam)
        sortIndexMinimap(sam2bam.out.bam)
        medaka(sortIndexMinimap.out.bamBai)

        getPhredTrimmed(nanoq.out.trimmedFastq.join(getReadNumberLength.out.numlenJSON))
        getMappedPercent(minimap2SE.out.sam.join(getPhredTrimmed.out.phredJSON))
        getBAMinfo(sortIndexMinimap.out.bamBai.join(getMappedPercent.out.mappedJSON))

        nextcladeDB()
        nextcladeSubtype(medaka.out.consensus, nextcladeDB.out.database)
        sierra(medaka.out.consensus.join(nextcladeSubtype.out.tsv))

        htmlRunReport(getBAMinfo.out.depthJSON.collect(), sierra.out.json.collect{it[0]})
        pdfReport(sierra.out.json, params.markdownFile)
}