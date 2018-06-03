import os

from text_handler import logger


#input: path string to current database of fastq files, provided by NGS process
#        path string to output all the results of MIXCR procedure
def mixcr_procedure(path_to_mixcr, fastq_path, outpath, chains, mmu):

    align_cmd, assemble_cmd, exportAlignments_cmd, exportClones_cmds = get_mixcr_cmds(fastq_path, outpath, mmu, chains)

    logger.info('Current wd is: ' + os.getcwd())
    logger.info('Changing wd to mixcr\'s dir')
    os.chdir(path_to_mixcr)
    logger.info('Current wd is: ' + os.getcwd())

    logger.info('Starting mixcr align procedure')
    cmd_exe(align_cmd)

    #execute assemble to clones command
    logger.info('Starting assemble command')
    cmd_exe(assemble_cmd)

    #execute export alignments command
    logger.info('Starting exportAlignments command')
    cmd_exe(exportAlignments_cmd)

    #execute export clones command for each chain separately
    for exportClones_cmd in exportClones_cmds:
        cmd_exe(exportClones_cmd)



'''assign variables to commands'''
def get_mixcr_cmds(fastq_path, outpath, MMU, chains):

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    logger.debug('fastq path: ' + fastq_path)
    logger.debug('os.path.join(fastq_path, \'R1.fastq\'): ' + os.path.join(fastq_path, 'R1.fastq'))

    fastq1 = fastq2 = ''
    for file_name in os.listdir(fastq_path):
        if 'R1.fastq' in file_name:
            fastq1 = os.path.join(fastq_path, file_name)
        elif 'R2.fastq' in file_name:
            fastq2 = os.path.join(fastq_path, file_name)
    logger.debug('fastq files paths are:\n{}\n{}'.format(fastq1, fastq2))

    if not os.path.exists(fastq1):
        logger.error('R1.fastq is missing...')
        raise OSError('R1.fastq does not exist...')

    if not os.path.exists(fastq2):
        logger.error('R2.fastq is missing...')
        raise OSError('R2.fastq does not exist...')

    vdjca_path = os.path.join(outpath, 'alignments.vdjca')
    clones_clns_path = os.path.join(outpath, 'clones.clns')

    align_cmd = ('java -Xmx4g -Xms3g -jar mixcr.jar align'                      #align command
                 ' -f'                                                          #overwrite output file if already exists
                 ' -s {}'                                                       #consider species (mouse/human)
                 ' --report ' + outpath + '/align_report.txt'                   #create report file
                 ' --library imgt'                                              #use IMGT local library as annotation reference
                 ' -a'                                                          #save reads' ids from fastq files
                 ' ' + fastq1 + ' ' + fastq2 + ''                               #input files- 2 X fastq files
                 ' ' + vdjca_path).format("mmu" if MMU else "hsa")
             
    assemble_cmd = ('java -Xmx4g -Xms3g -jar mixcr.jar assemble'                    #assemble command
                    ' -r ' + outpath + '/assemble_report.txt'                       #create report file
                    ' -f'                                                           #overwrite output file if already exists
                    ' -i ' + outpath + '/index_file'           #keep mapping between initial reads and final clones
                    ' -OseparateByC=true'                                           #separate by isotypes
                    ' -OcloneFactoryParameters.vParameters.featureToAlign=VRegion'  #align v region and not v transcript
                    ' -OassemblingFeatures=[CDR3]'                   #define sequence to create clones by
                    ' -OminimalClonalSequenceLength=6'                              #minimum number of nucleotides in clonal sequence
                    ' ' + vdjca_path  +      #input file - VDJCA from previous step
                    ' ' + clones_clns_path)                 #output file

    exportAlignments_cmd = ('java -Xmx4g -Xms3g -jar mixcr.jar exportAlignments'    #exportAlignments command   
                            ' -f'                                                   #overwrite output file if already exists
                            ' --preset-file aln_fields.txt'                         #export fields specified in aln_fields file
                            ' -cloneIdWithMappingType ' + outpath + '/index_file'  #indicate stase of each read
                            ' ' + vdjca_path  +   #input file- VDJCA from previous step
                            ' ' + outpath + '/alignments.txt')      #output file

    exportClones_cmds = []
    for chain in chains:
        exportClones_cmd = ('java -Xmx4g -Xms3g -jar mixcr.jar exportClones'            #exportClones command
                        ' -f'                                                           #overwrite output file if already exists
                        ' --chains ' + chain + ''
                        ' --preset-file assemble_fields.txt'                            #export fields specified in assemble_fields file
                        ' -readIds ' + outpath + '/index_file'
                        ' -o'                                                           #remove out-of-frame clones
                        ' -t'                                                           #remove stop codon clones
                        ' ' + clones_clns_path + ''                  #input file
                        ' ' + outpath + '/' + chain + '_clones.txt')
        exportClones_cmds.append(exportClones_cmd)

    return align_cmd, assemble_cmd, exportAlignments_cmd, exportClones_cmds


'''executing command in command prompt'''
def cmd_exe(cmd):
    logger.info('Executing:\n' + cmd)
    os.system(cmd)
        

                                