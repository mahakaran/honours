#Written by Mahakaran Sandhu, 2018-2019. 
#Honours Year, RSC, ANU. 

import operator
import numpy as np

def average(lst):
    return sum(lst) / len(lst)
    
def aa_PR_reader(string, alignment_length):
    """reads string output of PAML for a given site, returns dict of probability values for each AA. 
    #str --> dict"""
    i = 17+alignment_length+2+2
    j = 17+alignment_length+2+7
    aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    out_dict = {}
    for x in range(len(aa_list)):
     
        out_dict.update({aa_list[x]:float(string[i:j])})
        i = i+9
        j = j+9
        
    return out_dict       
    #tested, works successfully
    
def sequence_Pr(datalist, alignment_length):
    """reads PAML data for probabilities at a particular notes for the entire sequence, and returns a list of tuples, 
    the first item in the tuple is AA site, and the second item is the dictionary of probabilities as obtained by the
    function aa_Pr_reader. List --> list of tuples, where tuple = (str, dict)"""
    outlist = []
    for i in datalist: 
        site = i[2:10]
        
        aa_probabilities = aa_PR_reader(i, alignment_length )
        out_tuple = (site, aa_probabilities)
        outlist.append(out_tuple)
    return outlist
    
 def total_Pr(site_list):
    #takes the output of sequence_Pr to calculate 1) the most probable sequence and 2) the overall posterior probability
    #of the best possible sequence. List --> tuple(str, float)
    seq_list = []
    probabilities = []
    for i in site_list: 
        site = i[0]
        aa_probabilities = i[1]
        
        best_AA = max(aa_probabilities.iteritems(), key=operator.itemgetter(1))[0]
        best_AA_Pr = max(aa_probabilities.iteritems(), key=operator.itemgetter(1))[1]
        seq_list.append(best_AA)
        probabilities.append(best_AA_Pr)
    sequence_Pr = average(probabilities)
    best_sequence =  "".join(seq_list)
    return (best_sequence, sequence_Pr)

def position_frequency_matrix(site_list, node_name): 
    """takes the output of sequence_Pr and returns a PFM (numpy array). You can then print this array to a file for input
    #into a AA logo program. Outputs a file in TRANSFAC matrix format. """
    outfile = open(str(node_name)+'frequencies.transfac', 'w+')
    sites = range(len(site_list))
    aa_list = [('A',[]), ('R',[]), ('N',[]),
               ('D',[]),('C',[]), ('Q',[]), ('E',[]), ('G',[]), ('H',[]),('I',[]),
               ('L',[]), ('K',[]), ('M',[]), ('F',[]), ('P',[]), ('S', []), ('T', []),
               ('W',[]), ('Y',[]), ('V', [])]
 
    for i in site_list:
        aa_frequencies = i[1]
        for j in aa_list:
            freq = aa_frequencies.get(j[0])
            j[1].append(freq)
    
    freq_list = []
    freq_list.append(sites)
    for i in aa_list:
        freq_list.append(i[1])
        
    freq_array = np.array(freq_list)
    
    outfile.write('P0' + '\t A' +'\t R'+'\t N'+'\t D'+'\t C' + '\t Q' + '\t E' + '\t G' + '\t H' + '\t I' + '\t L'
                 + '\t K' + '\t M'+ '\t F' + '\t P' + '\t S'  + '\t T' + '\t W' + '\t Y' + '\t V' + '\n' )
    for i in range(len(freq_list[0])):
        outfile.write(str(i) + '\t' + str(freq_list[1][i]) + '\t' + str(freq_list[2][i]) + '\t' +str(freq_list[3][i])
                      +'\t'+str(freq_list[4][i])  + '\t'+str(freq_list[5][i])  + '\t'+str(freq_list[6][i])  + '\t'+
                      str(freq_list[7][i])  + '\t'+str(freq_list[8][i]) + '\t' +str(freq_list[9][i]) + '\t' +
                      str(freq_list[10][i]) + '\t'+str(freq_list[11][i])+ '\t' +str(freq_list[12][i]) + '\t' +
                      str(freq_list[13][i]) + '\t'+str(freq_list[14][i])+ '\t' +str(freq_list[15][i])  + '\t'+
                      str(freq_list[16][i]) + '\t' +str(freq_list[17][i])+ '\t'+str(freq_list[18][i])  + '\t'+
                      str(freq_list[19][i]) + '\t' +str(freq_list[20][i])+'\n')
        
    
    return freq_list
                
        
    
def read_rst(infile, node_list, seq_length):
    """reads the PAML output file rst, and pulls the correct lines out, i.e. the lines containing posterior probabilties
    when particular nodes are specified. Sequence length also needs to be specified. Writes a list of tuple, with 
    each tuple (node number, node data). """
    data = open(infile, 'r')
    datas = data.readlines()
    #if list element Prob dist at node x exists, append x number of following elements to a new list. 
    cat_data = []
    for i in range(len(datas)):
       
        
        node_data = []
        if 'Prob distribution at node' in datas[i]:
            if int(datas[i][26:29]) in node_list:
                node_data.append(datas[i+4:i+4+int(seq_length)])
                node = datas[i][26:29]
                node_data = node_data[0]
                nodes = (node, node_data)
                cat_data.append(nodes)
    return cat_data

def ASR_analyser(infile, node_list, seq_length, alignment_length):
    """for each node, run all the analysis functions defined above, and write correct output files. """
    ASR_data = read_rst(infile, node_list, seq_length)
    for i in ASR_data: 
        node = i[0]
        data = i[1]
        
        site_probabilities = sequence_Pr(data, alignment_length)
      
        best_sequence = total_Pr(site_probabilities)
        PFM = position_frequency_matrix(site_probabilities, node)
       
        outfile = open(str(node)+'_best_sequence', 'w+')
        outfile.write('>'+ str(node) + '\n')
        outfile.write(str(best_sequence[0] + '\n'))
        print node
        print best_sequence[1]
        

        
    
            
        
            
            
            
   
    
    
