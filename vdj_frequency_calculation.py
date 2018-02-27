import matplotlib     # Must be before importing matplotlib.pyplot or pylab! to Avoid the need of X-Dislay https://stackoverflow.com/questions/4706451/how-to-save-a-figure-remotely-with-pylab/4706614#4706614
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab! to Avoid the need of X-Dislay
import matplotlib.pyplot as plt
import pandas as pd
from collections import Counter
import re
from pylab import cm
from mpl_toolkits.mplot3d import Axes3D
import gzip


v_index = 2
d_index = 3
j_index = 4
start_index = 14
super_fam_indicator = [True, False]

def my_open(file_path, mode):
    if file_path.endswith('gzip'):
        return gzip.open(file_path, mode=mode)
    else:
        return open(file_path, mode=mode)

'''get VDJ families of current sequence'''       
def extract_family_assignment(line, chain, super_fam):
    
    #split description line to data pieces
    parsed_line = line.split("|")
    
    #extract super/sub family 
    v_fam = parsed_line[v_index][start_index:-1]
    if super_fam:
        v_fam = re.search('^[A-Za-z]+[0-9]+', v_fam).group()
    

    #extract d family in case of heavy chain
    if chain == "IGH":
        d_fam = parsed_line[d_index][start_index:-1]

        if d_fam != '' and super_fam: 
            d_fam = re.search('^[A-Za-z]+[0-9]+', d_fam).group()
    else:
        d_fam = ""

    j_fam = parsed_line[j_index][start_index:-1]
    if super_fam:
        j_fam = re.search('^[A-Za-z]+[0-9]+', j_fam).group()

    return v_fam, d_fam, j_fam
        

'''calculate VDJ frequency'''   
def get_vdj_combination_frequencies(comb_list):
 
    #get number of combinations (with duplications)
    sum_of_combinations = len(comb_list)

    #convert list with duplications to counter dictionary {"v d j" combination:count}
    vdj_comb_freq = Counter(comb_list)     
    
    #convert counts to frequencies
    for comb in vdj_comb_freq:
        vdj_comb_freq[comb] = round(( vdj_comb_freq[comb] / sum_of_combinations ) * 100 , 2)    

    return vdj_comb_freq

    
'''3D VDJ combination graph'''  
def plot_vdj_freq(vdj_freq_dict, sample, chain, super_fam, outpath):
    
    v_fam_set, d_fam_set, j_fam_set = set(), set(), set()
    v_cord_list, d_cord_list, j_cord_list, freq_list = list(), list(), list(), list()  
    
    #create list of all V, D, J families
    for vdj_str in vdj_freq_dict:
        vdj_comb = tuple(vdj_str.split())
        
        #VDJ combination => add families to set to 
        if len(vdj_comb) == 3:    
            v_fam_set.add(vdj_comb[0]), d_fam_set.add(vdj_comb[1]), j_fam_set.add(vdj_comb[2])
            
    #convert sets into ordered lists
    v_labels = list(sorted(v_fam_set))
    d_labels = list(sorted(d_fam_set))
    j_labels = list(sorted(j_fam_set))    

    #get coordinations of each vdj combination    
    for vdj_str in vdj_freq_dict:
        vdj_comb = tuple(vdj_str.split())
        
        #VDJ combination => get coordinates
        if len(vdj_comb) == 3:     
            v_cord_list.append(v_labels.index(vdj_comb[0]))
            d_cord_list.append(d_labels.index(vdj_comb[1]))
            j_cord_list.append(j_labels.index(vdj_comb[2]))
            freq_list.append(vdj_freq_dict[vdj_str])
            

    #plot the data            
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(111,projection='3d')
    cm.jet(freq_list)
    colmap = plt.get_cmap('jet', len(set(freq_list)))
    
    #size of combination depends on its frequency
    sizes = [freq_list[i]*100 for i in range(len(freq_list))]
    
    yg = ax.scatter(v_cord_list, d_cord_list, j_cord_list, c=freq_list, marker='o', s=sizes, cmap = colmap, alpha=0.5)

    #colorbar parameters
    cbar = fig.colorbar(yg)
    cbar.set_label("Frequency")
    
    #space axis labels from values
    ax.xaxis.labelpad=30
    ax.yaxis.labelpad=30    
    ax.zaxis.labelpad=30
    
    ax.set_xlabel('V family')
    ax.set_ylabel('D family')
    ax.set_zlabel('J family')
    
    plt.locator_params(axis='x',nbins=len(v_labels))
    plt.locator_params(axis='y',nbins=len(d_labels))
    plt.locator_params(axis='z',nbins=len(j_labels))
        
    ax.set_xlim([0,len(v_labels)])
    ax.set_ylim([0,len(d_labels)])
    ax.set_zlim([0,len(j_labels)])
    
    ax.set_xticklabels(v_labels)
    ax.set_yticklabels(d_labels)
    ax.set_zticklabels(j_labels)
    
    if super_fam:
        fig_path = outpath + "/" + sample + "/" + chain + "_VDJ_super_family.png"
    else:
        fig_path = outpath + "/" + sample + "/" + chain + "_VDJ.png"
        
    ax.grid(color='b')
    plt.title(sample + " " + chain + " VDJ combinations' frequencies")
    plt.savefig(fig_path)
    plt.close('all')
    
    
'''get combination of pairs (vd, vj, dv)'''    
def get_pair_combination_frequencies(vdj_freq_dict, all_fam_comb_set, first_fam, second_fam):
    
    pair_freq_dict = dict()            
    
    #iterate over vdj combinations
    for vdj in vdj_freq_dict:
        comb = tuple(vdj.split())
        
        #combinations contains "d"
        if len(comb) == 3:
            
            if first_fam == "v" and second_fam == "d":
                comb = (comb[0], comb[1])
            
            elif first_fam == "v" and second_fam == "j":
                comb = (comb[0], comb[2])
                
            elif first_fam == "d" and second_fam == "j":
                comb = (comb[1], comb[2])
        
        #combination does not contain "d"
        elif len(comb) == 2:
            
            if first_fam == "d" or second_fam == "d":
                continue
                
        #add chars_pair combination to frequency dictionary and save it in combinations' set
        if comb in pair_freq_dict: 
            pair_freq_dict [comb] += vdj_freq_dict[vdj]
        else:
            pair_freq_dict [comb] = vdj_freq_dict[vdj]
            all_fam_comb_set.add(comb)
    
    return pair_freq_dict 
    

'''get family assignments as list''' 
def get_families_list(all_fam_comb_set):
    
    #lists of all v/j families    
    a_families = set()
    b_families = set()
    
    for comb_tup in all_fam_comb_set:
        a_families.add(comb_tup[0])    
        b_families.add(comb_tup[1])    
    
    return a_families, b_families


'''get minimum and maximum frequency'''
def get_min_max_frequencies(all_smpl_chn_comb_freq, samples, chains):
    
    #initialize min/max values
    vmin = 101
    vmax = -1
    
    for sample in samples:
        for chain in chains:
            
            #get min/max frequency of sample and chain
            crnt_comb_dict = all_smpl_chn_comb_freq[(sample,chain)]
                                                    
            if len(crnt_comb_dict) == 0:
                continue
            
            crnt_min = min(crnt_comb_dict.values())
            crnt_max = max(crnt_comb_dict.values())
            
            #update new min/max if found
            if crnt_min < vmin:
                vmin = crnt_min
            if crnt_max > vmax:
                vmax = crnt_max
                
    return vmin, vmax


'''plot dataframe as square heatmap'''   
def plot_df(a_families, b_families, df, crnt_path, sample, chain, vmin, vmax, super_fam, pair):

    fig, ax = plt.subplots()
    
    #set axis positions and labels
    ax.set_xticks([i for i in range(len(a_families))])
    ax.set_yticks([i for i in range(len(b_families))])
    ax.set_xticklabels(sorted(a_families), rotation = 90)
    ax.set_yticklabels(sorted(b_families))
    
    #set plot parameters
    ax.tick_params(axis='both', labelsize = 180/max(len(a_families), len(b_families)))
    ax.grid(False)
    cmap_custom = plt.get_cmap('Reds')
    im = plt.imshow(df, interpolation='nearest', cmap = cmap_custom, vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(im)
    cbar.set_label("Frequency")
    plt.title(sample + " " + chain + " Frequency of " + pair + " combinations")
    plt.xlabel(pair[0] + " families")
    plt.ylabel(pair[1] + " families")
    
    #save image
    if super_fam:
        fig_path = crnt_path + "/" + chain + "_" + pair + "_super_family.png"
    else:
        fig_path = crnt_path + "/" + chain + "_" + pair + ".png"
        
    plt.savefig(fig_path)
    plt.close('all')

       
'''gather data for square heatmap'''  
def create_comb_heatmap(all_fam_comb_set, all_smpl_chn_comb_freq, super_fam, pair, samples, chains, outpath):
    
    #lists of chars_pair families
    a_families, b_families = get_families_list(all_fam_comb_set)        

    #initialize frequency matrix of vj combinations
    frequency_mat = [[0 for i in range(len(a_families))] for j in range(len(b_families))]

    #get min/max frequencies    
    vmin, vmax = get_min_max_frequencies(all_smpl_chn_comb_freq, samples, chains)
                        
    for sample in samples:
        for chain in chains:

            crnt_path = outpath + "/" + sample
            
            #fill frequency matrix according to frequency dictionary
            l=0
            for b_fam in sorted(b_families):
                k=0
                for a_fam in sorted(a_families):
                    if (a_fam, b_fam) in all_smpl_chn_comb_freq[(sample,chain)]:        
                        frequency_mat[l][k] = all_smpl_chn_comb_freq[(sample,chain)][(a_fam, b_fam)]
                    
                    k += 1
                l+=1
                        
            #convert maxtix to data frame
            df = pd.DataFrame(frequency_mat, index=b_families, columns=a_families)
            plot_df(a_families, b_families, df, crnt_path, sample, chain, vmin, vmax, super_fam, pair)
  

'''specify VDJ combination frequency in xls file'''
def extract_data_to_xls_file(vdj_freq_dict, sample, chain, sf, outpath):

    #set file path
    if sf:
        xls_path = outpath + "/" + sample + "/" + sample + "_" + chain + "_V(D)J_combination_frequency_super_family.xls"
    else:
        xls_path = outpath + "/" + sample + "/" + sample + "_" + chain + "_V(D)J_combination_frequency.xls"

    with open(xls_path, "w") as file:        

        file.write(sample + " " + chain + "V(D)J combination frequency\n")
        file.write("V family\tDfamily\tJfamily\tfrequency\n")
        
        for comb in sorted(vdj_freq_dict):
            
            comb_tup = comb.split()
            
            if len(comb_tup) == 3:
                for fam in comb_tup:
                    file.write(fam + "\t")
            elif len(comb_tup) == 2:
                file.write(comb_tup[0]+"\t\t"+comb_tup[1]+"\t")
                
            file.write(str(vdj_freq_dict[comb]))
                    
            file.write("\n")
    

'''VDJ combination analysis procedure'''
def get_genes_combination_freq(samples, chains, fasta_file, outpath):

    #repeat for super and sub families        
    for sf in super_fam_indicator:
        
        #dictionary of dictionaries. {(sample,chain):{(v,j):freq}}
        vd_all_smpl_chn_comb_freq = dict()
        vj_all_smpl_chn_comb_freq = dict()
        dj_all_smpl_chn_comb_freq = dict()
        
        #all combinations on vj familes
        vd_all_fam_comb_set = set()
        vj_all_fam_comb_set = set()
        dj_all_fam_comb_set = set()
        
        for sample in samples:
            for chain in chains:
                
                #list of all combinations
                comb_list = list()
                
                #initialize families
                v_fam = ""
                d_fam = ""
                j_fam = ""
                
                #iterate over sequences that their CDR3 appears more than once

                with my_open(fasta_file, "r") as f:
                    
                    res = -1
                    
                    for line in f:
                        
                        #get descriptve lines only
                        res += 1
                        if res%2 != 0:
                            continue
        
                        #make sure sequence belongs to current sample and chain
                        if "(" + sample + "," not in line:
                            continue
                        if "chain: " + chain not in line:
                            continue
                        
                        #get VDJ families and append combination to list
                        v_fam, d_fam, j_fam = extract_family_assignment(line, chain, sf)
                        comb = v_fam + " " + d_fam + " " + j_fam                        
                        comb_list.append(comb)
                        
                #frequency dictionary of all vdj combinations as strings
                vdj_freq_dict = get_vdj_combination_frequencies(comb_list)

                #create 3D plot of vdj combinations
                if chain == "IGH":
                    plot_vdj_freq(vdj_freq_dict, sample, chain, sf, outpath)
                    
                extract_data_to_xls_file(vdj_freq_dict, sample, chain, sf, outpath)
                    
                #frequency dictionary of all chars_pair combinations
                vj_freq_dict = get_pair_combination_frequencies(vdj_freq_dict, vj_all_fam_comb_set, "v","j")
                vd_freq_dict = get_pair_combination_frequencies(vdj_freq_dict, vd_all_fam_comb_set, "v","d")
                dj_freq_dict = get_pair_combination_frequencies(vdj_freq_dict, dj_all_fam_comb_set, "d","j")
                
                #add dictionary of chars_pair frequencies to sample and chain dictionary
                vj_all_smpl_chn_comb_freq[(sample,chain)] = vj_freq_dict
                vd_all_smpl_chn_comb_freq[(sample,chain)] = vd_freq_dict
                dj_all_smpl_chn_comb_freq[(sample,chain)] = dj_freq_dict

        #cretae heatmap of chars_pair combinations of all samples and chains
        create_comb_heatmap(vj_all_fam_comb_set, vj_all_smpl_chn_comb_freq, sf, "VJ", samples, chains, outpath)
        create_comb_heatmap(vd_all_fam_comb_set, vd_all_smpl_chn_comb_freq, sf, "VD", samples, chains, outpath)
        create_comb_heatmap(dj_all_fam_comb_set, dj_all_smpl_chn_comb_freq, sf, "DJ", samples, chains, outpath)
        
        
