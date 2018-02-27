import os
import re

cdr_index = 1
entry_id = 1


'''edit description line to specify only the sample that qualified the sequence'''
def edit_desc(sample, desc):
    
    parse_index = desc.rfind("|")

    #split desription line by last "|" to separate sample specification
    desc_start = desc[:parse_index+2]
    desc_end = desc[parse_index+2:]
    
    #number of times the sequence occur in the relevant sample
    seq_cnt = int( re.search("\("+sample+",(.*?)\)", desc_end).group(1) )
    
    #specify count for relevant sample
    edited_desc = desc_start + "(" + sample + ","
    
    return edited_desc, seq_cnt

'''provide dictionary of CDR3 counts and their sequences'''
def get_entries_dict(fasta_file, sample):

    with open(fasta_file, "r") as fas:
        
        #{CDR3:[CDR3 count, dictionary{(description,sequence):sequence count}}
        entries_dict = dict()
        
        #iterate over entries in the fasta file
        for desc_line in fas:
            
            desc_line = desc_line.strip()
            desc = desc_line[desc_line.index("chain:"):]
            desc, cdr_read_cnt = edit_desc(sample, desc)            
            seq = fas.readline().strip()
            cdr = desc_line.split("|")[cdr_index]

            if cdr in entries_dict:
                
                #update read count of CDR3
                entries_dict[cdr][0] += cdr_read_cnt

                #new entry with current CDR3 => add to set
                if (desc, seq) not in entries_dict[cdr][1]:
                    entries_dict[cdr][1][(desc, seq)] = cdr_read_cnt
                 
                #duplicate entry => keep one copy and update reads count
                else:
                   entries_dict[cdr][1][(desc, seq)] += cdr_read_cnt

            else:
                entries_dict[cdr] = [cdr_read_cnt, dict()]
                entries_dict[cdr][1][(desc, seq)] = cdr_read_cnt 
            
    return entries_dict 


'''write qualified sequecnes to final FASTA'''   
def write_qualified_sequences(final_fasta_path, all_seq_set, entries_dict, th):

    global entry_id
    
    with open(final_fasta_path , "a") as fas:
        
        for cdr in entries_dict:
            
            #CDR3 count of sequence above required threshold
            if entries_dict[cdr][0] >= th:
                
                #iterate over all sequences with current CDR3
                for entry in entries_dict[cdr][1]:
                    
                    #avoid sequence duplications in file
                    if entry[1] not in all_seq_set:
                        
                        #write description and sequence
                        fas.write("> " + str(entry_id) + " " + entry[0] + str(entries_dict[cdr][1][entry]) + ")\n")
                        fas.write(entry[1] + "\n")
                        entry_id += 1
                
                        #add sequence to set to avoid sequence duplications
                        all_seq_set.add(entry[1])
                        
 
'''write final FASTA by merging identical reads of different samples'''
def write_final_fasta(final_fasta_path_aux, final_fasta_path):
    
    all_seq_dict = dict()

    #iterate over all entries in auxiliary FASTA
    with open(final_fasta_path_aux, "r") as file:
        
        for desc_line in file:
            desc_line = desc_line.strip()
            seq = file.readline()
            
            #duplicate sequence in different samples => merge description lines 
            if seq in all_seq_dict:
                all_seq_dict[seq] += desc_line[desc_line.rfind("|")+1:]
                
            #sequence's first occurance
            else:
                all_seq_dict[seq] = desc_line
                
    #write merged entries to final FASTA file
    with open(final_fasta_path, "w") as file:
        for seq in all_seq_dict:
            file.write(all_seq_dict[seq] + "\n" + seq)

               
    
'''create FASTA of high confidence sequences'''
def create_final_fasta(samples, overlap_count_thresholds, non_overlap_count_thresholds, db_outpath):
    
    #path to fasta file with final sequences
    final_fasta_path = db_outpath + "/final_fasta.fasta"
    if os.path.exists(final_fasta_path):
        os.remove(final_fasta_path)
        
    #path to fasta file with final sequences - auxiliary
    final_fasta_path_aux = db_outpath + "/final_fasta_aux.fasta"
    if os.path.exists(final_fasta_path_aux):
        os.remove(final_fasta_path_aux)

        
    for i in range(len(samples)):

        #set of all written sequences, to avoid duplications
        all_seq_set = set()
        
        #original files' paths to extract sequences from
        sample_path = db_outpath + "/" + samples[i]
        overlap_path = sample_path + "/overlap_db.fasta"
        run1_path = sample_path + "/run1_db.fasta"
        run2_path = sample_path + "/run2_db.fasta"
        original_files = [overlap_path, run1_path, run2_path]
        
        for j in range(len(original_files)):

            if os.path.exists(original_files[j]):
                entries_dict = get_entries_dict(original_files[j], samples[i])

                #get required count threshold
                if j == 0:
                    th = overlap_count_thresholds[i]
                else:
                    th = non_overlap_count_thresholds[i]

                write_qualified_sequences(final_fasta_path_aux, all_seq_set, entries_dict, th)
    
    write_final_fasta(final_fasta_path_aux, final_fasta_path)


