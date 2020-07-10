import os
import subprocess
import shutil
import json
import ASAP_CONSTANTS as CONSTS
from text_handler import logger

def get_taxon_dict_from_json(lib_json, taxon):
    for taxon_dict in lib_json:
        if taxon_dict['taxonId'] == taxon:
            return taxon_dict

    raise ValueError('Cant really happen...')


def get_taxa_and_lib(lib_path):
    with open(lib_path) as f:
        lib_json = json.load(f)
    taxon2taxon_dict = {}
    for taxon_dict in lib_json:
        taxon = taxon_dict['taxonId']
        taxon2taxon_dict[taxon] = get_taxon_dict_from_json(lib_json, taxon)
    return taxon2taxon_dict, lib_json


def get_gene_names(alternative_taxon_dict):
    gene_names = set()
    for gene_dict in alternative_taxon_dict['genes']:
        gene_names.add(gene_dict['name'])
    return gene_names


def remove_colliding_entries_from_default_dict(alternative_taxon_dict, default_taxon_dict):
    alternative_gene_names = get_gene_names(alternative_taxon_dict)
    default_gene_names = set(get_gene_names(default_taxon_dict))
    for gene_name in alternative_gene_names:
        if gene_name in default_gene_names:
            gene_to_remove = -1
            for i, gene_dict in enumerate(default_taxon_dict['genes']):
                if gene_dict['name'] == gene_name:
                    gene_to_remove = i
                    break
            sequence_to_remove = -1
            for i, sequence_dict in enumerate(default_taxon_dict['sequenceFragments']):
                if sequence_dict['uri'].endswith(gene_name):
                    sequence_to_remove = i
                    break
            if gene_to_remove > -1:
                if sequence_to_remove == -1:
                    logger.error(f'Did not find sequence for {gene_name}')
                default_taxon_dict['genes'].pop(gene_to_remove)
                default_taxon_dict['sequenceFragments'].pop(sequence_to_remove)


def remove_collisions_from_default_lib(local_default_lib_path, alternative_lib_path):
    alternative_taxon2taxon_dict, alternative_lib_json = get_taxa_and_lib(alternative_lib_path)
    default_taxon2taxon_dict, default_lib_json = get_taxa_and_lib(local_default_lib_path)
    for taxon in alternative_taxon2taxon_dict:
        if taxon in default_taxon2taxon_dict:
            remove_colliding_entries_from_default_dict(alternative_taxon2taxon_dict[taxon], default_taxon2taxon_dict[taxon])

    result = [value for value in default_taxon2taxon_dict.values()]
    # save manipulated local default lib
    with open(local_default_lib_path, 'w') as f:
        json.dump(result, f, indent=4)


# remove_collisions_from_default_lib('/Users/Oren/Desktop/imgt.201822-5.sv4 2.json',
#                                    '/Users/Oren/Desktop/alternative_lib 2.json')


def generate_lib_for_mixcr(wd, default_lib_path, alternative_lib_path, output_html_path, remote_run):
    try:
        # avoid changing the original default lib!!
        local_default_lib_path = f'{wd}/manipulated_default_lib.json'
        shutil.copyfile(default_lib_path, local_default_lib_path)
        remove_collisions_from_default_lib(local_default_lib_path, alternative_lib_path)

        merged_lib_path = f'{wd}/merged_lib.json'
        execute_command(f'repseqio merge -f {local_default_lib_path} {alternative_lib_path} {merged_lib_path}')
        logger.info(f'Merge succeeded! Using merged lib at: {merged_lib_path}')
        return merged_lib_path
    except:
        # merging failed. will use the default lib
        logger.info(f'Merging failed! Using ORIGINAL default lib at: {default_lib_path}')
        if remote_run:
            with open(output_html_path) as f:
                html_txt = f.read()
            logger.info(f"html_txt.find('+ alternative_lib'): {html_txt.find('+ alternative_lib')}")
            html_txt = html_txt.replace('+ alternative_lib', '(merging failed)').replace('refr-lib', 'red-text')
            with open(output_html_path, 'w') as f:
                f.write(html_txt)
        return default_lib_path

#input: path string to current database of fastq files, provided by NGS process
#        path string to output all the results of MIXCR procedure
def mixcr_procedure(fastq_path, outpath, mmu, lib_path, remote_run, error_path):

    align_cmd, assemble_cmd, exportAlignments_cmd, exportClones_cmds = get_mixcr_cmds(lib_path, fastq_path, outpath, mmu, remote_run, error_path)

    logger.info('Current wd is: ' + os.getcwd())
    #logger.info('Changing wd to mixcr\'s dir')
    #os.chdir(path_to_mixcr)
    #logger.info('Current wd is: ' + os.getcwd())

    #for debugging:
    library_cmd = 'mixcr -v'.split()  #which library paths mixcr uses
    logger.info(f'Germline library paths of mixcr are:\n{subprocess.check_output(library_cmd).decode()}')

    logger.info('Starting mixcr align procedure')
    execute_command(align_cmd)

    #execute assemble to clones command
    logger.info('Starting assemble command')
    execute_command(assemble_cmd)

    #execute export alignments command
    logger.info('Starting exportAlignments command')
    execute_command(exportAlignments_cmd)

    #execute export clones command for each chain separately
    #TODO: why it is needeD? Answer: it is NOT.
    #for exportClones_cmd in exportClones_cmds:
    #    execute_command(exportClones_cmd)



'''assign variables to commands'''
def get_mixcr_cmds(lib_path, fastq_path, outpath, MMU, remote_run, error_path):

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    logger.debug(f'fastq path: {fastq_path}')
    logger.debug(f'os.path.join(fastq_path, "R1.fastq"): {os.path.join(fastq_path, "R1.fastq")}')

    fastq1 = fastq2 = ''
    for file_name in os.listdir(fastq_path):
        if 'fastq' in file_name:
            if 'R1' in file_name:
                fastq1 = os.path.join(fastq_path, file_name)
            elif 'R2' in file_name:
                fastq2 = os.path.join(fastq_path, file_name)
    logger.info(f'fastq files paths are:\n{fastq1}\n{fastq2}')

    if not os.path.exists(fastq1):
        logger.error('R1.fastq is missing...')
        raise OSError('R1.fastq does not exist...')

    if not os.path.exists(fastq2):
        logger.error('R2.fastq is missing...')
        raise OSError('R2.fastq does not exist...')

    verify_fastq_files_format(error_path, fastq1, fastq2)

    vdjca_path = os.path.join(outpath, 'alignments.vdjca')
    clones_clns_path = os.path.join(outpath, 'clones.clns')

    align_cmd = ('mixcr align'                      #align command
                 ' -f'                                                          #overwrite output file if already exists
                 f' -s {"mouse" if MMU else "human"}'                                                       #consider species (mouse/human)
                 ' -c IGH,IGL,IGK'                                              #immunological chain gene(s) to align
                 #f' --report {outpath}/align_report.txt'                   #create report file
                 f' --library {lib_path.split(".json")[0]}'  # mixcr requires lib name without json suffix!!
                 ' -a'                                                          #save reads' ids from fastq files
                 ' --threads 4'                                                          #number of threads
                 #' --verbose'
                 f' {fastq1} {fastq2}'                               #input files- 2 X fastq files
                 f' {vdjca_path}')
             
    assemble_cmd = ('mixcr assemble'                    #assemble command
                    ' -r ' + outpath + '/assemble_report.txt'                       #create report file
                    ' -f'                                                           #overwrite output file if already exists
                    f' -i {outpath}/index_file'                                   #keep mapping between initial reads and final clones
                    ' -OseparateByC=true'                                           #separate by isotypes
                    ' --threads 4'                                                          #number of threads
                    #' -OcloneFactoryParameters.vParameters.featureToAlign=VRegion' #align v region and not v transcript
                    #' -OassemblingFeatures=[CDR3]'                   #define sequence to create clones by
                    #' -OminimalClonalSequenceLength=6'                             #minimum number of nucleotides in clonal sequence
                    f' {vdjca_path}'                                             #input file - VDJCA from previous step
                    f' {clones_clns_path}')                                         #output file

    exportAlignments_cmd = ('mixcr exportAlignments'    #exportAlignments command   
                            ' -f'                                        #overwrite output file if already exists
                            f' --preset-file {CONSTS.ASAP_EXEC+"/" if remote_run else ""}aln_fields.txt'         #export fields specified in aln_fields file
                            f' -cloneIdWithMappingType {outpath}/index_file'     #indicate stase of each read
                            f' {vdjca_path}'                                   #input file- VDJCA from previous step
                            f' {outpath}/alignments.txt')                           #output file

    exportClones_cmds = []
    # for chain in chains:
    #     exportClones_cmd = ('mixcr exportClones'                #exportClones command
    #                     ' -f'                                   #overwrite output file if already exists
    #                     f' --chains {chain}'
    #                     f' --preset-file {CONSTS.ASAP_EXEC+"/" if remote_run else ""}assemble_fields.txt'    #export fields specified in assemble_fields file
    #                     f' -readIds {outpath}/index_file'
    #                     ' -o'                                   #remove out-of-frame clones
    #                     ' -t'                                   #remove stop codon clones
    #                     f' {clones_clns_path}'                  #input file
    #                     f' {outpath}/{chain}_clones.txt')
    #     exportClones_cmds.append(exportClones_cmd)

    return align_cmd, assemble_cmd, exportAlignments_cmd, exportClones_cmds


def verify_fastq_files_format(error_path, fastq1, fastq2):

    rep_num = os.path.split(os.path.split(fastq1)[0])[-1][-1] #/bioseq/data/results/asap/154832296135203243128690777655/reads/run1/R1.fastq

    num_lines_fastq1 = sum(1 for line in open(fastq1) if line.rstrip() != '')
    logger.info(f'{num_lines_fastq1} lines in {fastq1}')

    if num_lines_fastq1 % 4 != 0:
        err_msg = f'Illegal fastq file format: number of lines in {os.path.split(fastq1)[-1]} of rep {rep_num} is not a multiple of 4. One or more records are faulty.'
        logger.error(err_msg)
        with open(error_path, 'w') as error_path_f:
            error_path_f.write(err_msg)
        raise ValueError(err_msg)

    num_lines_fastq2 = sum(1 for line in open(fastq2) if line.rstrip() != '')
    logger.info(f'{num_lines_fastq2} lines in {fastq2}')

    if num_lines_fastq2 % 4 != 0:
        err_msg = f'Illegal fastq file format: number of lines in {os.path.split(fastq2)[-1]} of rep {rep_num} is not a multiple of 4. One or more records are faulty.'
        logger.error(err_msg)
        with open(error_path, 'w') as error_path_f:
            error_path_f.write(err_msg)
        raise ValueError(err_msg)

    if num_lines_fastq1 != num_lines_fastq2:
        err_msg = f'Illegal fastq files format: {os.path.split(fastq1)[-1]} and {os.path.split(fastq2)[-1]} of rep {rep_num} contain different number of lines ({num_lines_fastq1} and {num_lines_fastq2}, respectively).'
        logger.error(err_msg)
        with open(error_path, 'w') as error_path_f:
            error_path_f.write(err_msg)
        raise ValueError(err_msg)


'''executing command in command prompt'''
def execute_command(cmd):
    logger.info(f'Executing mixcr:\n{cmd}')
    logger.info(subprocess.check_output(cmd, shell=True).decode('unicode_escape'))
        

                                