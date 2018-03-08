import os
import subprocess

from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio.Align.Generic import Alignment
from Bio.Alphabet import Gapped, IUPAC

from mixcr_procedure import mixcr_procedure
from parse_alignments import parse_alignment_file
from directory_creator import create_dir
from runs_aggregator import join_runs_analyses
from plots_generator import plot_barplot
from text_handler import write_dict_to_file, string_similarity
from msa_parser import remove_sparse_columns
import logging

from weblogo_generator import generate_weblogo

logger = logging.getLogger('main')

def analyze_samples(gp):

    for i in range(1, gp.number_of_runs+1):

        run = 'run' + str(i)

        mixcr_output_path, parsed_mixcr_output_path, assignments_path = create_sub_working_directories(gp, run)

        try:
            done_path = os.path.join(gp.output_path, 'done_1_' + run + '_' + 'mixcr_procedure.txt')
            if os.path.exists(done_path):
                logger.info('Skipping mixcr_procedure, output files already exist...')
            else:
                logger.info('Starting mixcr_procedure of {}...'.format(run))
                fastq_path = os.path.join(gp.working_dir, 'reads', run)
                mixcr_procedure(gp.path_to_mixcr, fastq_path, mixcr_output_path, gp.sample, gp.chains, gp.MMU)
                with open(done_path, 'w') as f:
                    pass
        except:
            logger.error('Error in MiXCR procedure')
            raise

        try:
            done_path = os.path.join(gp.output_path, 'done_2_' + run + '_' + 'parse_alignment_file.txt')
            if os.path.exists(done_path):
                logger.info('Skipping parse_alignment_file, output files already exist...')
            else:
                logger.info('Applying parse alignments procedure of {}...'.format(run))
                parse_alignment_file(mixcr_output_path, parsed_mixcr_output_path, gp.sequence_annotation_file_suffix, gp.mutation_count_file_suffix, gp.len_threshold, gp.sequencing_quality_threshold)
                with open(done_path, 'w') as f:
                    pass
        except:
            logger.error('Error in parse alignments procedure')
            raise


    try:
        logger.info('Joining runs...')
        '''
        joint_run_is_needed = True #TODO: should it be always OR when gp.run_number > 1?
        if joint_run_is_needed:
        '''
        mixcr_output_path, parsed_mixcr_output_path, assignments_path = create_sub_working_directories(gp, 'joint')
        gp.parsed_mixcr_output_paths.append(parsed_mixcr_output_path)
        gp.assignments_paths.append(assignments_path)

        done_path = os.path.join(gp.output_path, 'done_3_' + 'join_runs_analysis.txt')
        if os.path.exists(done_path):
            logger.info('Skipping join_runs_analysis, output files already exist...')
        else:
            join_runs_analyses(gp.number_of_runs, gp.run_output_paths, parsed_mixcr_output_path, gp.chains,
                               gp.sequence_annotation_file_suffix, gp.mutation_count_file_suffix)
            with open(done_path, 'w') as f:
                pass
    except:
        logger.error('Error in join_runs_analysis')
        raise

    try:
        done_path = os.path.join(gp.output_path, 'done_4_' + 'parse_annotation_files.txt')
        if os.path.exists(done_path):
            logger.info('Skipping parse_annotation_files, output files already exist...')
        else:
            parse_annotation_files(gp)
            with open(done_path, 'w') as f:
                pass
    except:
        logger.error('Error in parse_annotation_files')
        raise

    try:
        done_path = os.path.join(gp.output_path, 'done_5_' + 'plot_assignments.txt')
        if os.path.exists(done_path):
            logger.info('Skipping plot_assignments, output files already exist...')
        else:
            plot_assignments(gp)
            with open(done_path, 'w') as f:
                pass
    except:
        logger.error('Error in plot_assignments')
        raise

    try:
        done_path = os.path.join(gp.output_path, 'done_6_' + 'plot_mutation_counts.txt')
        if os.path.exists(done_path):
            logger.info('Skipping plot_mutation_counts, output files already exist...')
        else:
            plot_mutation_counts(gp)
            with open(done_path, 'w') as f:
                pass
    except:
        logger.error('Error in plot_mutation_counts')
        raise

    try:
        done_path = os.path.join(gp.output_path, 'done_7_' + 'analyze_cdr3.txt')
        if os.path.exists(done_path):
            logger.info('Skipping analyze_cdr3, output files already exist...')
        else:
            cdr3_analysis_dir = os.path.join(gp.output_path, 'joint', 'cdr3_analysis')
            create_dir(cdr3_analysis_dir)
            analyze_cdr3(gp, parsed_mixcr_output_path, cdr3_analysis_dir)
            with open(done_path, 'w') as f:
                pass
    except:
        logger.error('Error in analyze_cdr3')
        raise


def plot_mutation_counts(gp):
    for parsed_mixcr_output_path in gp.parsed_mixcr_output_paths:
        # if gp.debug_run not in parsed_mixcr_output_path:
        #     logger.info('!!DEBUG!! SKIPPING parse_chain_annotation_file for', parsed_mixcr_output_path)
        #     continue
        for chain in gp.chains:
            mutation_counts_frequency_path = parsed_mixcr_output_path + '/' + chain + gp.mutation_count_file_suffix
            logger.info('Plotting mutation_counts_frequency_file: ' + mutation_counts_frequency_path)
            if os.path.exists(mutation_counts_frequency_path):
                plot_barplot(mutation_counts_frequency_path, gp.raw_data_file_suffix, key_type=int, value_type=int, fontsize=6, rotation=70, ylim=[0,30])
            else:
                logger.info('No such file: {}. Skipping this plot...'.format(mutation_counts_frequency_path))


def plot_assignments(gp):
    for assignments_path in gp.assignments_paths:
        if gp.debug_run not in assignments_path:
            logger.info('!!DEBUG!! SKIPPING parse_chain_annotation_file for ' + assignments_path)
            continue
        for raw_data_file in os.listdir(assignments_path):
            if not raw_data_file.endswith('_counts.' + gp.raw_data_file_suffix):
                logger.debug('Skipping ' + raw_data_file)
                continue
            assignment_file = os.path.join(assignments_path, raw_data_file)

            logger.info('Plotting assignment_file: ' + assignment_file)
            if assignment_file.endswith('cdr3_len_counts.' + gp.raw_data_file_suffix):
                # TODO: maybe this file should be handled separately
                plot_barplot(assignment_file, gp.raw_data_file_suffix, key_type=int)
            else:
                if 'VDJ' in assignment_file:
                    plot_barplot(assignment_file, gp.raw_data_file_suffix, rotation=70, fontsize=2, ylim=[0, 15])
                elif 'VD' in assignment_file or 'DJ' in assignment_file or 'VJ' in assignment_file:
                    plot_barplot(assignment_file, gp.raw_data_file_suffix, rotation=70, fontsize=5, ylim=[0, 25])
                else:
                    plot_barplot(assignment_file, gp.raw_data_file_suffix)


def parse_annotation_files(gp):
    for i in range(gp.number_of_runs + 1):  # +1 for the joint analysis (if needed)
        if i == gp.number_of_runs:
            run = 'joint run'
        else:
            run = 'run' + str(i + 1)
        if gp.debug_run not in run:
            logger.info('!!DEBUG!! SKIPPING parse_chain_annotation_file for', run)
            continue
        logger.info('Starting parse chain annotation file of ' + run)
        parse_sequence_annotation_file(gp.parsed_mixcr_output_paths[i], gp.assignments_paths[i],
                                       gp.sequence_annotation_file_suffix, gp.raw_data_file_suffix)


def analyze_cdr3(gp, parsed_mixcr_output_path, cdr3_analysis_dir):
    k = gp.top_cdr3_clones
    for chain in gp.chains:
        cdr3_to_counts = {}
        cdr3_to_aa_reads = {}
        most_k_common_cdr3_to_entry = {}

        annotations_path = os.path.join(parsed_mixcr_output_path, chain + gp.sequence_annotation_file_suffix)
        with open(annotations_path) as f:
            for line in f:
                aa_read, chain, cdr3, v_type, d_type, j_type, dna_read, isotype, read_frequency = line.split('\t')
                cdr3_to_counts[cdr3] = cdr3_to_counts.get(cdr3, 0) + int(read_frequency)
                cdr3_to_aa_reads[cdr3] = cdr3_to_aa_reads.get(cdr3, []) + [aa_read]
        most_k_common_cdr3 = sorted(cdr3_to_counts, key=cdr3_to_counts.get, reverse=True)[:k]
        most_k_common_cdr3_to_counts = {cdr3: cdr3_to_counts[cdr3] for cdr3 in most_k_common_cdr3}
        most_k_common_cdr3_to_aa_reads = {cdr3: cdr3_to_aa_reads[cdr3] for cdr3 in most_k_common_cdr3}
        if most_k_common_cdr3 != []:
            for i in range(k):
                cdr3 = most_k_common_cdr3[i]
                multiple_sequences = most_k_common_cdr3_to_aa_reads[cdr3]
                cluster_prefix = os.path.join(cdr3_analysis_dir, 'cluster_' + str(i))
                ms_path = cluster_prefix + '_ms.fasta'
                msa_path = cluster_prefix + '_msa.aln'

                try:
                    aln, msa = align_sequences(multiple_sequences, ms_path, msa_path, i)
                except:
                    print(multiple_sequences)
                    raise
                summary_align = AlignInfo.SummaryInfo(aln)

                # find consensus sequence by majority rule
                consensus = summary_align.dumb_consensus(threshold=0, ambiguous='#')
                consensus = str(consensus)

                #find most similar sequence
                aa_most_similar_to_consesnsus, similarity_rate = find_most_similar_sequence(msa, consensus)

                dna_most_similar_to_consesnsus = find_correspnding_dna(aa_most_similar_to_consesnsus, annotations_path)

                #add entry to dict
                most_k_common_cdr3_to_entry[cdr3] = [cdr3, str(most_k_common_cdr3_to_counts[cdr3]), consensus, aa_most_similar_to_consesnsus, dna_most_similar_to_consesnsus, str(similarity_rate)]

                #generate web logo
                weblogo_file_path = cluster_prefix + '_weblogo.pdf'
                try:
                    generate_weblogo(msa_path, weblogo_file_path)
                except Exception as e:
                    logger.warning('Couldn\'t generate weblogo {} due to the following reason: {}'.format(msa_path, e.args))

            #write cdr3_to_entry to file
            cdr3_annotations_path = os.path.join(cdr3_analysis_dir, chain + gp.cdr3_annotation_file_suffix)
            with open(cdr3_annotations_path, 'w') as f:
                for cdr3 in most_k_common_cdr3:
                    f.write('\t'.join(most_k_common_cdr3_to_entry[cdr3]) + '\n')
            #write_dict_to_file(cdr3_annotations_path, most_k_common_cdr3_to_entry, sort_by=most_k_common_cdr3_to_counts.get, reverse=True)


def find_most_similar_sequence(msa, consensus):
    max_similarity = 0
    max_seq = ''
    for seq in msa:
        similarity = round(string_similarity(seq, consensus), 2)
        if similarity > max_similarity:
            max_similarity = similarity
            max_seq = seq
    return max_seq, max_similarity


def find_correspnding_dna(aa_most_similar_to_consesnsus, annotations_path):
    # assuming that the dna sequence appears at the 6th position
    '''
    b'QVQLVQSGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSGISWNSGSIGYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTALYYCAKDLIVAVPAAAKVSAFDIWGQGTMVTVSS\tIGH\tCAKDLIVAVPAAAKVSAFDIW\tIGHV3\tIGHD6\tIGHJ3\tCCCAGGTCCAGCTGGTGCAGTCTGGGGGAGGCTTGGTACAGCCTGGCAGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTGATGATTATGCCATGCACTGGGTCCGGCAAGCTCCAGGGAAGGGCCTGGAGTGGGTCTCAGGTATTAGTTGGAATAGTGGTAGCATAGGCTATGCGGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCCCTGTATCTGCAAATGAACAGTCTGAGAGCTGAGGACACGGCCTTGTATTACTGTGCAAAAGATCTTATTGTAGCAGTACCAGCTGCCGCAAAAGTGAGTGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAGGGAGTGCATCCGCCCCAACCCTCTCTCTCTCTCCTCCGGG\t6\n'
    '''
    raw_dna_most_similar_to_consesnsus = subprocess.check_output(['grep', aa_most_similar_to_consesnsus, annotations_path]).split()[6]
    dna_most_similar_to_consesnsus = str(raw_dna_most_similar_to_consesnsus)[2:-1] #without ^b' an $'
    logger.info('dna_most_similar_to_consesnsus is {}'.format(dna_most_similar_to_consesnsus))
    return str(dna_most_similar_to_consesnsus)


def parse_sequence_annotation_file(parsed_mixcr_output_path, assignments_path, sequence_annotation_file_suffix, raw_data_file_suffix, allowed_chain_types = ['IGH', 'IGK', 'IGL']):

    for chain in allowed_chain_types:

        chain_data_path = os.path.join(parsed_mixcr_output_path, chain + sequence_annotation_file_suffix)

        if not os.path.exists(chain_data_path):
            logger.info('Skipping: ' + chain_data_path + ' (file does not exists)')
            continue

        v_types = {}
        d_types = {}
        j_types = {}
        vd_types = {}
        vj_types = {}
        dj_types = {}
        vdj_types = {}
        cdr3_len = {}
        cdr3_frequency_counter = {}

        with open(chain_data_path) as f:

            for line in f:
                aa_read, chain, cdr3, v_type, d_type, j_type, dna_read, isotype, read_frequency = line.rstrip().split('\t')

                v_types[v_type] = v_types.get(v_type, 0) + 1
                d_types[d_type] = d_types.get(d_type, 0) + 1
                j_types[j_type] = j_types.get(j_type, 0) + 1
                vd_types[v_type+' '+d_type] = vd_types.get(v_type+' '+d_type, 0) + 1
                vj_types[v_type+' '+j_type] = vj_types.get(v_type+' '+j_type, 0) + 1
                dj_types[d_type+' '+j_type] = dj_types.get(d_type+' '+j_type, 0) + 1
                vdj_types[v_type+' '+d_type+' '+j_type] = vdj_types.get(v_type+' '+d_type+' '+j_type, 0) + 1
                cdr3_len[len(cdr3)] = cdr3_len.get(len(cdr3), 0) + 1
                cdr3_frequency_counter[cdr3] = cdr3_frequency_counter.get(cdr3, 0) + int(read_frequency)

            v_path = assignments_path + '/' + chain + '_V_counts.' + raw_data_file_suffix
            write_dict_to_file(v_path, v_types)

            d_path = assignments_path + '/' + chain + '_D_counts.' + raw_data_file_suffix
            write_dict_to_file(d_path, d_types)

            j_path = assignments_path + '/' + chain + '_J_counts.' + raw_data_file_suffix
            write_dict_to_file(j_path, j_types)

            vd_path = assignments_path + '/' + chain + '_VD_counts.' + raw_data_file_suffix
            write_dict_to_file(vd_path, vd_types)

            vj_path = assignments_path + '/' + chain + '_VJ_counts.' + raw_data_file_suffix
            write_dict_to_file(vj_path, vj_types)

            dj_path = assignments_path + '/' + chain + '_DJ_counts.' + raw_data_file_suffix
            write_dict_to_file(dj_path, dj_types)

            vdj_path = assignments_path + '/' + chain + '_VDJ_counts.' + raw_data_file_suffix
            write_dict_to_file(vdj_path, vdj_types)

            cdr3_len_path = assignments_path + '/' + chain + '_cdr3_len_counts.' + raw_data_file_suffix
            write_dict_to_file(cdr3_len_path, cdr3_len)

            write_dict_to_file(parsed_mixcr_output_path + '/cdr3_frequency_counter.' + raw_data_file_suffix, cdr3_frequency_counter, sort_by=cdr3_frequency_counter.get, reverse=True)


def create_sub_working_directories(global_params, run):

    run_output_path = os.path.join(global_params.output_path, run)
    global_params.run_output_paths.append(run_output_path)
    create_dir(run_output_path)

    mixcr_output_path = ''
    if run != 'joint':
        mixcr_output_path = os.path.join(run_output_path, 'mixcr')
        global_params.mixcr_output_paths.append(mixcr_output_path)
        create_dir(mixcr_output_path)

    parsed_mixcr_output_path = os.path.join(run_output_path, 'parsed_mixcr_output')
    create_dir(parsed_mixcr_output_path)
    global_params.parsed_mixcr_output_paths.append(parsed_mixcr_output_path)

    assignments_path = os.path.join(run_output_path, 'vdj_assignments')
    create_dir(assignments_path)
    global_params.assignments_paths.append(assignments_path)

    return mixcr_output_path, parsed_mixcr_output_path, assignments_path


def align_sequences(multiple_sequences, ms_path, msa_path, cluster_number):
    # create alignment object from a list of multiple sequences and write it to a file

    aln = Alignment(Gapped(IUPAC.protein, '-'))

    # write ms to file
    with open(ms_path, 'w') as f:
        f.write('\n'.join('>'+str(i) + '\n' + multiple_sequences[i] for i in range(len(multiple_sequences))))

    if len(multiple_sequences) > 1:
        # align file
        logger.info('running mafft:')
        mafft_cmd = 'mafft ' + ms_path + ' > ' + msa_path
        logger.info(mafft_cmd)
        os.system(mafft_cmd)
        # remove columns with more than 20% gaps
        remove_sparse_columns(msa_path, msa_path, 0.2)
    else:
        logger.info('Only 1 sequence in {}. Skipping mafft and copying ms file as is.'.format(str(cluster_number)))
        with open(msa_path, 'w') as f:
            f.write('\n'.join('>' + str(i) + '\n' + multiple_sequences[i] for i in range(len(multiple_sequences))))

    for record in SeqIO.parse(msa_path, 'fasta'):
        aln.add_sequence(record.id, str(record.seq))

    logger.debug('Alignment is:')
    msa_str = ''
    for record in aln._records:
        msa_str += record.seq._data + '\n'
    logger.debug(msa_str)

    return aln, msa_str
