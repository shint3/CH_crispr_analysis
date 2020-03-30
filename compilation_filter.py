import csv
import sys

#type into command line on Biowulf:
#python compilation_filter.py ref_seq_with_IUPAC allele_table_1 ... allele_table_n
#
#ref_seq_with_IUPAC = reference sequence with IUPAC
#                     (to account for heterozygous SNPs unique to a specific macaque)
#allele_table_n = alleles frequency table from CRISPResso


#list of canonical nucleotides
#(for determining which nucleotides in the reference seuqence are noncanonical
#and need to be replaced with IUPAC)
nt_lst=['A','T','G','C']

#dictionary for saving noncanonical nucleotides and their positions in the reference sequence
iupac_nts = {}

#variable for saving the command line argument number of the alleles frequency table
#that is currently being used
allele_table_num = 2

#dictionary of read numbers at different timepoints for each unique read sequence
read_num_dict = {}
#dictionary of read percentages at different timepoints for each unique read sequence
read_percent_dict = {}

#variable for saving the position in the list of read numbers/percentages
#within the dictionary of read numbers/percentages that is currently being updated
timepoint_num = 0

#list of zeroes that is used to construct
#list of read numbers at different timepoints for a unique read sequence
read_num_zeroes_lst = []
#list of zeroes that is used to construct
#list of read percentages at different timepoints for a unique read sequence
read_percent_zeroes_lst = []

#list of column titles (names of alleles frequency tables)
col_titles = []

#threshold value for the filter
#(for a read sequence to be kept, at least one timepoint's read percentage after editing
#needs to differ by more than the threshold value from the read percentage befor editing)
filter_threshold = 0.1
#variable used to count how many times
#the read percentage of a sequence at a timepoint after editing differs by more than
#the threshold value from the read percentage of the sequence before editing
similar_count = 0

#for each command line argument
for arg in sys.argv:

    #if the 2nd arguemnt (if the reference sequence with IUPAC)
    if sys.argv.index(arg) == 1:

        #for each nucleotide in reference sequence that is noncanonical (not A, T, G, or C),
        #save their position as a key and IUPAC notation as the corresponding value
        #in the iupac_nts dictionary
        iupac_ref = arg
        for ref_nt_num,ref_nt in enumerate(iupac_ref):
            if ref_nt not in nt_lst:
                iupac_nts[ref_nt_num] = ref_nt

    #if the 3rd or later argument (if an alleles frequency table)
    if sys.argv.index(arg) > 1:

        #open alleles frequency table
        with open(sys.argv[allele_table_num],"U") as allele_table:
            reader = csv.reader(allele_table, delimiter ='\t')

            #for each row in the alleles frequency table
            for row_num,allele_table_row in enumerate(reader):

                #skip 1st row with column titles
                if row_num == 0:
                    continue

                #for all later rows with data (sequence, read number, read percentage),
                #change the nucleotide(s) in the sequence to IUPAC notation
                #if those nucleotide(s) are IUPAC notation in the reference sequence
                #(with this info having been previously saved in iupac_nts dictionary)
                else:
                    row_seq = list(allele_table_row[0])                       
                    for seq_nt_num,seq_nt in enumerate(row_seq):
                        if seq_nt_num in iupac_nts:
                            row_seq[seq_nt_num] = iupac_nts[seq_nt_num]
                        allele_table_row[0] = ''.join(row_seq)

                #if the read sequence (after IUPAC) is the same as the reference sequence,
                #or if CRISPResso has identified an insertion or deletion in the sequence
                if allele_table_row[0] == sys.argv[1] or allele_table_row[3] =='False':

                    #if the read sequence is not already in the dictionary of read numbers,   
                    if allele_table_row[0] not in read_num_dict:

                        #make a list of zeroes equal in length to the number of alleles frequency tables    
                        for num in range(2,len(sys.argv)):
                            read_num_zeroes_lst.append(0)

                        #add a new key to the dictionary of read numbers equal to the read sequence
                        #and change the corresponding value to be equal to the list of zeroes
                        #then update the list to include the read number at the sequence's timepoint
                        read_num_dict[allele_table_row[0]] = read_num_zeroes_lst
                        read_num_dict[allele_table_row[0]][timepoint_num] = int(allele_table_row[8])


                        #make a list of zeroes equal in length to the number of alleles frequency tables
                        for num in range(2,len(sys.argv)):
                            read_percent_zeroes_lst.append(0)

                        #add a new key to the dictionary of read percents equal to the read sequence
                        #and change the corresponding value to be equal to the list of zeroes
                        #then update the list to include the read percent at the sequence's timepoint   
                        read_percent_dict[allele_table_row[0]] = read_percent_zeroes_lst
                        read_percent_dict[allele_table_row[0]][timepoint_num] = float(allele_table_row[9])

                        #empty lists of zeroes
                        read_num_zeroes_lst = []
                        read_percent_zeroes_lst = []

                    #if the read sequence is already in the dictionary of read numbers,
                    else:

                        #update the value of the read sequence key to include the read number at the new timepoint
                        read_num_dict[allele_table_row[0]][timepoint_num]= read_num_dict[allele_table_row[0]][timepoint_num] + int(allele_table_row[8])

                        #update the value of the read sequence key to include the read percentage at the new timepoint
                        read_percent_dict[allele_table_row[0]][timepoint_num]= read_percent_dict[allele_table_row[0]][timepoint_num] + float(allele_table_row[9])
                        
            #move on to working with next allele frequency table inputted in command line                              
            allele_table_num += 1
            #input data into the next column of the tables of read numbers/percentages
            timepoint_num += 1


#construct list of names of the alleles frequency tables
#(using a method independent of the file extension of the alleles frequency tables)
for num,arg in enumerate(sys.argv):
    if num > 1:
        col_title = arg.split('.')
        col_title = col_title[:-1]
        col_titles.append('.'.join(col_title))

#print out the list of alleles frequency table names as the column titles of new data file
print '\t'.join(col_titles)

#for each key in the dictionary of read percents
for read_seq in read_percent_dict:

    #compare read percentage of sequence before editing with
    #the read percentage of each timepoint after editing,
    #and if the percentage of a timepoint after editing differs
    #by more than the threshold value,
    #then increase the counter variable (similar_count) by 1
    for num,percent in enumerate(read_percent_dict[read_seq]):
        if num == 0:
            before_edit_percent = percent
        if num > 0:
            if abs(float(percent - before_edit_percent)) > filter_threshold:
                similar_count += 1

    #if the counter variable (similar_count) is greater than 0
    #(if at least one timepoint's percentage after editing differs
    #by more than the threshold value from the percentage before editing),
    #then print out the key's corresponding value in the dictionary of read numbers
    if similar_count > 0:
        print read_seq,'\t', '\t'.join(map(str,read_num_dict[read_seq]))

    #reset counter variable (similar_count) to 0
    similar_count=0
