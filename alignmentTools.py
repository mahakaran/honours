import Bio
from Bio import AlignIO
def remove_seqs(alignment,alnformat, removal_ids, outname): 

    full_align = AlignIO.read(alignment, alnformat)
    removal_list = open(removal_ids, 'r')
    removal_list = removal_list.readlines()
    removal_list = [i.strip() for i in removal_list]
    
    purged_ids = [record.id for record in full_align if record.id not in removal_list]
    purged_seqs = [record.seq for record in full_align if record.id not in removal_list]
    
    outfile = open(outname, 'w+')
    for i in range(len(purged_ids)): 
        outfile.write('>' + purged_ids[i] + '\n')
        outfile.write(str(purged_seqs[i])+'\n')
