import os

from mixcr_procedure import mixcr_procedure
from parse_alignments import parse_alignment_file
from directory_creator import create_dir
from runs_aggregator import join_runs_analyses
from plots_generator import plot_barplot, generate_mutations_boxplots
from text_handler import write_dict_to_file, read_table_to_dict
from cdr3_analyzer import analyze_cdr3

import logging

logger = logging.getLogger('main')

def analyze_samples(gp):

    for i in range(1, gp.number_of_runs+1):

        run = 'run' + str(i)

        mixcr_output_path, parsed_mixcr_output_path, assignments_path, cdr3_analysis_path = create_sub_working_directories(gp, run)

        try:
            done_path = os.path.join(gp.output_path, 'done_1_' + run + '_' + 'mixcr_procedure.txt')
            if os.path.exists(done_path):
                logger.info('Skipping mixcr_procedure, output files already exist...')
            else:
                logger.info('Starting mixcr_procedure of {}...'.format(run))
                fastq_path = os.path.join(gp.working_dir, 'reads', run)
                mixcr_procedure(gp.path_to_mixcr, fastq_path, mixcr_output_path, gp.chains, gp.MMU)
                with open(done_path, 'w') as f:
                    pass
        except Exception as e:
            logger.error('Error in MiXCR procedure: {}'.format(str(e)))
            raise

        try:
            done_path = os.path.join(gp.output_path, 'done_2_' + run + '_' + 'parse_alignment_file.txt')
            if os.path.exists(done_path):
                logger.info('Skipping parse_alignment_file, output files already exist...')
            else:
                logger.info('Applying parse alignments procedure of {}...'.format(run))
                parse_alignment_file(mixcr_output_path, parsed_mixcr_output_path, gp.sequence_annotation_file_suffix, gp.mutations_file_suffix, gp.len_threshold, gp.sequencing_quality_threshold)
                with open(done_path, 'w') as f:
                    pass
        except Exception as e:
            logger.error('Error in parse alignments procedure: {}'.format(str(e)))
            raise e

    remove_irrelevant_chains(gp)

    if gp.chains == []:
        logger.info('No relevant chains to analyze...')
        raise ValueError('No relevant chains to analyze...')

    try:
        logger.info('Joining runs...')
        # TODO: uncomment this block once final_fasta is implemented for each file
        '''
        joint_run_is_needed = gp.run_number > 1 
        if joint_run_is_needed:
        '''
        joint_mixcr_output_path, joint_parsed_mixcr_output_path, joint_assignments_path, joint_cdr3_analysis_path = create_sub_working_directories(gp, 'joint')

        done_path = os.path.join(gp.output_path, 'done_3_' + 'join_runs_analysis.txt')
        if os.path.exists(done_path):
            logger.info('Skipping join_runs_analysis, output files already exist...')
        else:
            logger.info('Joining runs analysis...')
            join_runs_analyses(gp.number_of_runs, gp.run_output_paths, joint_parsed_mixcr_output_path, gp.chains,
                               gp.sequence_annotation_file_suffix, gp.mutations_file_suffix)
            with open(done_path, 'w') as f:
                pass
    except Exception as e:
       logger.error('Error in join_runs_analyses: {}'.format(str(e)))
       raise

    try:
        done_path = os.path.join(gp.output_path, 'done_4_' + 'parse_annotation_files.txt')
        if os.path.exists(done_path):
            logger.info('Skipping parse_annotation_files, output files already exist...')
        else:
            logger.info('Parsing annotations files...')
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
            logger.info('Plotting V(D)J assignments...')
            plot_assignments(gp)
            with open(done_path, 'w') as f:
                pass
    except Exception as e:
       logger.error('Error in plot_assignments: {}'.format(str(e)))
       raise

    try:
        done_path = os.path.join(gp.output_path, 'done_6_' + 'plot_mutation_analyses.txt')
        if os.path.exists(done_path):
            logger.info('Skipping plot_mutation_analyses, output files already exist...')
        else:
            logger.info('Plotting mutation analyses...')
            plot_mutation_analyses(gp)
            with open(done_path, 'w') as f:
                pass
    except Exception as e:
       logger.error('Error in plot_mutation_counts: {}'.format(str(e)))
       raise

    try:
        done_path = os.path.join(gp.output_path, 'done_7_' + 'analyze_cdr3.txt')
        if os.path.exists(done_path):
            logger.info('Skipping analyze_cdr3, output files already exist...')
        else:
            logger.info('Analyzing CDR3s...')
            analyze_cdr3(gp)
            with open(done_path, 'w') as f:
                pass
    except Exception as e:
       logger.error('Error in analyze_cdr3: {}'.format(str(e)))
       raise


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

    cdr3_analysis_path = os.path.join(run_output_path, 'cdr3_analysis')
    create_dir(cdr3_analysis_path)
    global_params.cdr3_analysis_paths.append(cdr3_analysis_path)

    return mixcr_output_path, parsed_mixcr_output_path, assignments_path, cdr3_analysis_path


def remove_irrelevant_chains(gp):
    chains_involvement = dict.fromkeys(gp.chains, 0)
    for chain in gp.chains:
        for i in range(1, gp.number_of_runs + 1):
            parsed_mixcr_output_path = os.path.join(gp.output_path, 'run' + str(i), 'parsed_mixcr_output')
            chain_annotation_file = chain + gp.sequence_annotation_file_suffix
            if chain_annotation_file in os.listdir(parsed_mixcr_output_path):
                chains_involvement[chain] += 1
                break
    for chain in chains_involvement:
        if chains_involvement[chain] == 0:
            # no annotations were extracted for this chain. no point to further analyze it
            logger.error('No annotations were extracted for {}. Removing it from chains list...'.format(chain))
            gp.chains.remove(chain)


def parse_annotation_files(gp):
    for i in range(gp.number_of_runs + 1):  # +1 for the joint analysis (if needed)
        if i == gp.number_of_runs:
            run = 'joint run'
        else:
            run = 'run' + str(i + 1)
        logger.info('Starting parse chain annotation file of ' + run)
        parse_sequence_annotation_file(gp.parsed_mixcr_output_paths[i], gp.assignments_paths[i], gp.cdr3_analysis_paths[i],
                                       gp.sequence_annotation_file_suffix, gp.raw_data_file_suffix, skip_rows=1)


def parse_sequence_annotation_file(parsed_mixcr_output_path, assignments_path, cdr3_analysis_path, sequence_annotation_file_suffix, raw_data_file_suffix, allowed_chain_types = ['IGH', 'IGK', 'IGL'], skip_rows=0):
    for chain in allowed_chain_types:

        chain_data_path = os.path.join(parsed_mixcr_output_path, chain + sequence_annotation_file_suffix)

        if not os.path.exists(chain_data_path):
            logger.info('Skipping: ' + chain_data_path + ' (file does not exists)')
            continue

        v_types = {} #:{str:int}
        d_types = {} #:{str:int}
        j_types = {} #:{str:int}
        vd_types = {} #:{str:int}
        vj_types = {} #:{str:int}
        dj_types = {} #:{str:int}
        vdj_types = {} #:{str:int}
        cdr3_len = {} #:{int:int}
        cdr3_frequency_counter = {} #:{str:int}

        with open(chain_data_path) as f:
            for i in range(skip_rows):
                f.readline()
            for line in f:
                chain, isotype, dna_read, aa_read, cdr3, v_type, d_type, j_type, read_frequency = line.rstrip().split('\t')

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
            write_dict_to_file(v_path, v_types, header='\t'.join(['V family subgroup', 'counts']))

            d_path = assignments_path + '/' + chain + '_D_counts.' + raw_data_file_suffix
            write_dict_to_file(d_path, d_types, header='\t'.join(['D family subgroup', 'counts']))

            j_path = assignments_path + '/' + chain + '_J_counts.' + raw_data_file_suffix
            write_dict_to_file(j_path, j_types, header='\t'.join(['J family subgroup', 'counts']))

            vd_path = assignments_path + '/' + chain + '_VD_counts.' + raw_data_file_suffix
            write_dict_to_file(vd_path, vd_types, header='\t'.join(['VD combinations', 'counts']))

            vj_path = assignments_path + '/' + chain + '_VJ_counts.' + raw_data_file_suffix
            write_dict_to_file(vj_path, vj_types, header='\t'.join(['VJ combinations', 'counts']))

            dj_path = assignments_path + '/' + chain + '_DJ_counts.' + raw_data_file_suffix
            write_dict_to_file(dj_path, dj_types, header='\t'.join(['DJ combinations', 'counts']))

            vdj_path = assignments_path + '/' + chain + '_VDJ_counts.' + raw_data_file_suffix
            write_dict_to_file(vdj_path, vdj_types, header='\t'.join(['VDJ combinations', 'counts']))

            cdr3_len_path = cdr3_analysis_path + '/' + chain + '_cdr3_len_counts.' + raw_data_file_suffix
            write_dict_to_file(cdr3_len_path, cdr3_len, header='\t'.join(['length_of_CDR3', 'counts']))

            write_dict_to_file(cdr3_analysis_path + '/cdr3_frequency_counter.' + raw_data_file_suffix, cdr3_frequency_counter, sort_by=cdr3_frequency_counter.get, reverse=True, header='\t'.join(['CDR3', 'counts']))


def plot_assignments(gp):
    for assignments_path in gp.assignments_paths:
        for raw_data_file in os.listdir(assignments_path):
            if not raw_data_file.endswith(gp.raw_data_file_suffix):
                continue #avoid reading png files when re-running
            assignment_file = os.path.join(assignments_path, raw_data_file)
            plot_path = assignment_file.replace(gp.raw_data_file_suffix, 'png')
            logger.info('Plotting assignment_file: ' + assignment_file)
            d = read_table_to_dict(assignment_file, value_type=int, skip_rows=1)

            if 'VDJ' in assignment_file:
                plot_barplot(d, plot_path, x_label='\nV(D)J subgroups combinations')
                # plot_barplot(assignment_file, gp.raw_data_file_suffix)
            elif 'VD' in assignment_file:
                plot_barplot(d, plot_path, x_label='\nVD subgroups combinations')
            elif 'DJ' in assignment_file:
                plot_barplot(d, plot_path, x_label='\nDJ subgroups combinations')
            elif 'VJ' in assignment_file:
                plot_barplot(d, plot_path, x_label='\nVJ subgroups combinations')
            elif 'V' in assignment_file:
                plot_barplot(d, plot_path, x_label='\nV subgroups')
            elif 'D' in assignment_file:
                plot_barplot(d, plot_path, x_label='\nD subgroups')
            elif 'J' in assignment_file:
                plot_barplot(d, plot_path, x_label='\nJ subgroups')
            else:
                logger.error('Skipping non-assignment file: ' + raw_data_file)


def plot_mutation_analyses(gp):
    for parsed_mixcr_output_path in gp.parsed_mixcr_output_paths:
        for chain in gp.chains:

            mutations_path = parsed_mixcr_output_path + '/' + chain + gp.mutations_file_suffix
            if os.path.exists(mutations_path):
                dna_to_Ka_Ks_dict, dna_to_mutation_counts = parse_mutations_raw_data(mutations_path)
                logger.info('Plotting mutations_file: {}'.format(mutations_path))

                mutation_counts_to_frequency = get_values_frequency(dna_to_mutation_counts)
                mutation_counts_to_frequency_plot = mutations_path.replace(gp.raw_data_file_suffix, '1.png')
                plot_barplot(mutation_counts_to_frequency, mutation_counts_to_frequency_plot, x_label='\nNumber of basepair mutations')

                Ka_Ks_box_plot = mutations_path.replace(gp.raw_data_file_suffix, '2.png')
                generate_mutations_boxplots(dna_to_Ka_Ks_dict, Ka_Ks_box_plot)
            else:
                logger.info('Skipping plot for {} (no such file...)'.format(mutations_path))

            # mutation_counts_path = parsed_mixcr_output_path + '/' + chain + gp.mutations_file_suffix
            # if os.path.exists(mutation_counts_path):
            #     logger.info('Plotting mutation_counts_file: {}'.format(mutation_counts_path))
            #     #raw_data = read_table_to_dict(mutation_counts_path, key_type=int, value_type=int)
            #     mutation_counts_plot = mutation_counts_path.replace('txt', 'png')
            #     plot_barplot(mutation_counts_path, gp.raw_data_file_suffix, key_type=int, value_type=float, fontsize=6, ylim=[0,30])
            # else:
            #     logger.info('Skipping plot for {} (no such file...)'.format(mutation_counts_path))


            '''#old mutation plot before ka_ks
            mutation_counts_frequency_path = parsed_mixcr_output_path + '/' + chain + gp.mutations_file_suffix
            logger.info('Plotting mutation_counts_frequency_file: ' + mutation_counts_frequency_path)
            if os.path.exists(mutation_counts_frequency_path):
                plot_barplot(mutation_counts_frequency_path, gp.raw_data_file_suffix, key_type=int, value_type=float, fontsize=6, rotation=70, ylim=[0,30])
            else:
                logger.info('No such file: {}. Skipping this plot...'.format(mutation_counts_frequency_path))
            '''


def parse_mutations_raw_data(Ka_Ks_path):
    Ka_Ks_data = {} #:{str:[float,float]}
    mutations_data = {} #:{str:int}
    raw_data = read_table_to_dict(Ka_Ks_path, value_type=list, skip_rows=1)
    for key in raw_data:
        values = raw_data[key]
        Ka_Ks_data[key] = values[:2]
        mutations_data[key] = int(values[2])
    return Ka_Ks_data, mutations_data


def get_values_frequency(d):
    result = {}
    for value in d.values():
        result[value] = result.get(value, 0) + 1
    return result