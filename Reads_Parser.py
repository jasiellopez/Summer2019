import docx
import collections
import operator

def main():
#was supposed to be script which returned most common motifs that occurred within a given DNA sequence
#never fully developed it as we eventually lost the need for it
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

main()
