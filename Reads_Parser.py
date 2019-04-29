import docx
import collections
import operator

def main():

    files = input("Please input the files you wish to search, with a space between each file path")
    files_list = []
    files_list = files.split()
    sequences_dict = {}
    #dictionary where key is sequence name as input in the file name and value is sequence (sequence name based on FASTA format)
    files_dict = {}
    #dictionary where the key is the file and value is the sequences_dict corresponding to that file

    newfiles_list = []
    for filename in files_list:
        #/Users/jasiel/Downloads/Molgulid.mesp.alignments.edited_highlighted_motifs.docx
        #new_filename = filename + "_highlighted_motifs.docx"
        #new_file = docx.Document(new_filename)
        #newfiles_list.append(new_file)
        file = open(filename, "r")
        sequence_as_list = file.readlines()
        one_line = ""
        i = 0
        while i < len(sequence_as_list):
            if len(sequence_as_list[i]) > 0 and sequence_as_list[i][0] != '>':
                one_line += sequence_as_list[i].strip("\n")
                if len(sequence_as_list[i]) < line_length:
                    sequences_dict[sequences_name] = one_line
                i += 1
            elif sequence_as_list[i][0] == '>':
                #get sequence length
                sequences_name = sequence_as_list[i]
                split_sequence_name = sequences_name.split()
                #get line_length
                line_length = int(split_sequence_name[1].strip('\n'))
                i += 2
        files_dict[filename] = sequences_dict

    motifs_dict1 = {}
    start_site_dict = {}
    for filename in files_dict:
        motifs_dict1[filename] = {}
        start_site_dict[filename] = {}
        for sequence in files_dict[filename]:
            motifs_dict1[filename][sequence] = {}
            start_site_dict[filename][sequence] = 0

    for filename in files_dict:
        for name in files_dict[filename]:
            #name in this context is the sequence name
            sequence1 = files_dict[filename][name]
            print(sequence1)
            remainder = len(sequence1) % 4
            atgFound = False
            start_site = 0
            i = len(sequence1) - 1
            while i in range(len(sequence1)):
                if atgFound == True:
                    if (i - 4) >= 0:
                        potential_motif1 = sequence1[i - 4:i]
                        complement = GenerateComplement(potential_motif1)
                        CompSeq = False
                        inSeq = False
                        if potential_motif1 in motifs_dict1[filename][name] or complement in motifs_dict1[filename][name]:
                            if potential_motif1 in motifs_dict1[filename][name]:
                                if filename in motifs_dict1:
                                    if name in motifs_dict1[filename]:
                                        motifs_dict1[filename][name][potential_motif1].append(i)
                            if complement in motifs_dict1[filename][name]:
                                if filename in motifs_dict1:
                                    if name in motifs_dict1[filename]:
                                        motifs_dict1[filename][name][complement].append(i)
                        else:
                            motifs_dict1[filename][name][potential_motif1] = [i]
                    i -= 1
                else:
                    potential_site = sequence1[(i - 2):(i + 1)]
                    if (potential_site == "ATG"):
                        start_site = i - 2
                        start_site_dict[filename][name] = start_site
                        print(start_site)
                        atgFound = True
                    if atgFound == False:
                        i -= 1
                    elif atgFound == True:
                        i = start_site

    for filename in motifs_dict1:
        print("In file " + filename)
        for name in motifs_dict1[filename]:
            print("\tIn sequence " + name)
            num_motif_occurences = {}
            most_common_motifs = {}
            for motif in motifs_dict1[filename][name]:
                num_motif_occurences[motif] = len(motifs_dict1[filename][name][motif])
            for j in range(5):
                key = max(num_motif_occurences.items(), key=operator.itemgetter(1))[0]
                num_motif_occurences.pop(key)
                most_common_motifs[key] = motifs_dict1[filename][name][key]
            for motif in most_common_motifs:
                print("\t\tThe motif " + motif + " occurred at the following bp upstream of the start site:")
                indexes = ""
                for index in motifs_dict1[filename][name][motif]:
                    indexes += str(index) + "; "
                print("\t\t\t" + indexes)


def GenerateComplement(sequence):

    complement = ""
    i = len(sequence) - 1
    while i < len(sequence) and i >= 0:
        if sequence[i] == "A":
            complement += "T"
        elif sequence[i] == "G":
            complement += "C"
        elif sequence[i] == "C":
            complement += "G"
        elif sequence[i] == "T":
            complement += "A"
        i -= 1

    return complement

#find ATG site for protein coding sequence, and then work upstream scanning 4bp at a time while only moving one index at a time
main()

"ASK WHY MESP CIS-REG SEQUENCE (IN WORD DOC) SHOWS UP AS EXON ON ANISEED"



'''paragraphs = []
            #TODO: DO OTHER STUFF
            #inserts a -1 when a new alignment sequence is starteds
            for motif in motifs_found1:
                motifs_dict1[motif].append(-1)
            for motif in motifs_found2:
                motifs_dict2[motif].append(-1)
            motifs_found1 = []
            motifs_found2 = []
        for key in motifs_dict1:

    correction = 0
    for i in range(0, (len(sequence1)//60)):
        print((len(sequence1)//60))
        print("IN IF STATEMENT")
        print(i)
        print(len(sequences_dict[key][0]))
        corrected_i = (i*60) + correction
        print(corrected_i)
        pre_new_line = len(sequences_dict[key][0])
        print(sequences_dict[key][0][:50])
        print(sequences_dict[key][0][:54])
        print(sequences_dict[key][0][:58])
        print(sequences_dict[key][0][:62])
        sequences_dict[key][0] = sequences_dict[key][0][:corrected_i] + "\n" + sequences_dict[key][0][corrected_i:]
        print(len(sequences_dict[key][0]))
        print(sequences_dict[key][0][:corrected_i])

        print(sequences_dict[key][1])
        sequences_dict[key][1] = sequences_dict[key][1][:corrected_i] + "\n" + sequences_dict[key][1][corrected_i:]
        print("\n")
        print(sequences_dict[key][1])

        correction = len(sequences_dict[key][0]) - pre_new_line
'''
#add two lines at a time in a run (current line and alignment)

#find distance between 2 motifs



#consult andrew
#find overrpresentation close to an identified ETS site

#ask joe whether we shouldve been done with specs we ran












'''import docx

def main():

    motifs = input("Please input motifs you would like to check for, with a space in between each motif")
    motifs_dict = {}
    motifs_list = motifs.split()
    for motif in motifs_list:
        motifs_dict[motif] = []

    files = input("Please input the files you wish to search, with a space between each file path")
    files_list = []
    files_list = files.split()
    sequences_dict = {}

    for filename in files_list:
        file = open(filename, "r")
        sequence = file.readlines()
        one_line = ""
        alignment = ""
        i = 0
        while i < len(sequence):
            if len(sequence[i]) > 0 and sequence[i][0] != '>':
                one_line += sequence[i].strip("\n")
                alignment += sequence[i + 1].strip("\n")
                if len(sequence[i]) < 60:
                    sequences_dict[sequences_name] = [one_line, alignment]
                i += 3
            if i < len(sequence):
                if sequence[i][0] == '>':
                    sequences_name = sequence[i]
                    i += 2
    print(sequences_dict)

    newfiles_list = []
    for i in len(sequences_dict):
        for
        #/Users/jasiel/Downloads/Molgulid.mesp.alignments.edited_highlighted_motifs.docx
        new_filename = files_list[i] + "_highlighted_motifs.docx"
        new_file = docx.Document(new_filename)
        remainder = 0
        curr_line = 0
        #denotes whether previous line was divisible by 4
        previous_line = False
        for line in file:
            if ((len(line) + remainder) % 4 != 0):
                remainder = len(line) % 4
                for j in range(len(line)):
                    if j <= len(line) - remainder - 1:

                        potential_motif = line[j:j + 4]
                        #potential_motif =
                        if potential_motif in motifs_dict:
                            found_motif_index = j;
                            motifs_dict[potential_motif].append(found_motif_index)
                            #highlighttext
                            #increment occurence amount
                    if j > len(line) - remainder - 1:
                        sequence_remainder = line[j:]
                previous_line = True
            else:
                for j in range(len(line)):
                    potential_motif = line[j:j + 4]
                    if potential_motif in motifs_dict:
                        found_motif_index = j;
                        motifs_dict[potential_motif].append(found_motif_index)
                        #highlighttext
                        #increment occurence amount
                previous_line = False
            curr_line += 1

    print (motifs_dict)
#find distance between 2 motifs
'''
