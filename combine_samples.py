import os
import pickle
import re
from text_handler import read_table_to_dict
from global_params import sample


'''provide dictionary of all sequences and their counts in each sample'''
def combine_samples_count_dictionaries(outpath, unique_sequences_cnt):

    #initialize dictionary {sequence:[1,{sample_number:seq count in sample}]}
    count_seq_all_samples = dict()
        

    #load from file dictionary of sequences and their count
    crnt_dict = read_table_to_dict(os.path.join(outpath, 'sequences_frequency_counter.' + gp.raw_data_file_suffix + ''))

    #iterate over unique sequences of sample
    for seq in crnt_dict:

        #sequence belongs to more than one sample
        if seq in count_seq_all_samples:

            sample_count_dict = count_seq_all_samples[seq][1]   #get value of sequence, dictionary of sample ids and counts
            sample_count_dict[sample] = crnt_dict[seq]       #add sample id and its count

        #add sequence to dictionary
        else:
            count_seq_all_samples[seq] = [1,dict({sample: crnt_dict[seq]})]

    #WHY THE HELL DOES THIS NECESSARY??
    #compress dictionary when done with it
    #pickle.dump(crnt_dict, open(file_name, "wb"))

    return count_seq_all_samples 
       
      
'''create 3 FASTA files (all sequences, CDR3 repetitions, sequences repetitions)'''      
def combine_all_subsets(parsed_mixcr_output_path, all_samples_seq_count_dict):
    
    random_id = 0
    cdr3_set = set()
    
    all_file_non_single = open(parsed_mixcr_output_path + "/../all_subset_unique_sequences_non_single.fasta", "w")
    all_file_cdr3_rep = open(parsed_mixcr_output_path + "/../all_subset_unique_sequences_cdr3_rep.fasta", "w")
    
    #open file for sequences of all sub sets
    with open(parsed_mixcr_output_path + "/../all_subset_unique_sequences.fasta", "w") as all_file:

        for file_name in os.listdir(parsed_mixcr_output_path):

            if file_name in ['IGH.fasta', 'IGK.fasta', 'IGL.fasta']:

                with open(os.path.join(parsed_mixcr_output_path, file_name)) as f:

                    #iterate over entries in current file
                    for descriptive_line in f:

                        #AA sequence
                        seq_line = f.readline().strip()

                        #find sequence in dictionary
                        if seq_line in all_samples_seq_count_dict:

                            #sequence count information in all sub sets
                            sample_info_dict = all_samples_seq_count_dict[seq_line]

                            #sequence has not been seen yet => add it to file
                            if sample_info_dict[0] == 1:

                                #indicate that sequence has been seen => don't write it again
                                all_samples_seq_count_dict[seq_line][0] = 0

                                #extract sequence description
                                start_desc_index = descriptive_line.find("chain")

                                description = descriptive_line[start_desc_index:].strip()

                                #prepare string of all sub sets that hold current sequence with its count in the sub set
                                counts_dict = sample_info_dict[1]
                                count_str = ""
                                for sample, count in counts_dict.items():
                                    count_str += "(" + sample + "," + str(count) + ") "

                                #add sequence to full file
                                random_id += 1
                                all_file.write("> " + str(random_id) + " " + description + " | " + count_str + "\n")
                                all_file.write(seq_line + "\n")

                                #sequence has more than one appearance => add sequence to non single file
                                samples_cnt = len(counts_dict)
                                if samples_cnt > 1 or (samples_cnt == 1 and list(counts_dict.items())[0][1] != 1):
                                    all_file_non_single.write("> " + str(random_id) + " " + description + " | " + count_str + "\n")
                                    all_file_non_single.write(seq_line + "\n")

                                else:
                                    #find CDR3 of current sequence
                                    cdr3 = re.search('CDR3: (.+?) |', description)
                                    if cdr3:
                                        cdr3 = cdr3.group(1)

                                        #cdr3 appears in more than one sequence => add it to the file
                                        if cdr3 in cdr3_set:
                                            all_file_cdr3_rep.write("> " + str(random_id) + " " + description + " | " + count_str + "\n")
                                            all_file_cdr3_rep.write(seq_line + "\n")

                                        #cdr3 first appearance => add it to the set
                                        else:
                                            cdr3_set.add(cdr3)
                        else:
                            logger.info(50*'@' + 'seq line not in all_samples_seq_count_dict!!!!:\n', seq_line)



    all_file_non_single.close()
    all_file_cdr3_rep.close()
                                    
                 