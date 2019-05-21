import sys
import operator
import collections

def main():

    inputs = ParseArgs(sys.argv)
    motifs_input = inputs[0]
    files_list = inputs[1]
    window = inputs[2]
    motif_size = inputs[3]

    motifs_dict1 = {}
    motif_possibilities = {}
    motifs_generated = [] #list of all motifs generated from input motifs with n in them or variations (need to implement)
    start_site_dict = {}
    for motif in motifs_input:
        motif_possibilities[motif] = GenerateComplements(motif)
        motifs_dict1[motif] = {}

    sequences_dict = {} #dictionary where key is alignment and value is list of alignment
    files_dict = {} #dictionary where the key is the file and vaule is its alignments

    newfiles_list = []
    for filename in files_list:
        #/Users/jasiel/Desktop/Molgulid.mesp.alignments.edited_highlighted_motifs.docx
        #new_filename = filename + "_highlighted_motifs.docx"
        #new_file = docx.Document()
        #new_file.save(new_filename)
        #newfiles_list.append(new_file)
        file = open(filename, "r")
        sequence_as_list = file.readlines()
        '''
        num_times = 0
        p = new_file.add_paragraph()
        list_styles = new_file.styles
        highlighted = list_styles.add_style('highlighted', WD_STYLE_TYPE.CHARACTER)
        highlighted.base_style = list_styles['Default Paragraph Font']
        highlighted.font.highlight_color = WD_COLOR_INDEX.YELLOW
        string_sequence = ''.join(sequence)
        for motif in motifs_list:
            curr_motif = string_sequence.split(motif)
            for i in range(len(curr_motif) - 1):
                if num_times == 0:
                    p.add_run(curr_motif[i])
                    if len(curr_motif) > 1:
                        motif_run = p.add_run(text=motif, style=highlighted).font
                        motif_run.highlight_color = WD_COLOR_INDEX.YELLOW
        new_file.add_paragraph(p.text)
        new_file.save(new_filename)
        '''
        alignment = ""
        i = 0
        while i < len(sequence_as_list):
            if len(sequence_as_list[i]) > 0 and sequence_as_list[i][0] != '>':
                one_line += sequence_as_list[i].strip("\n")
                i += 1
            elif sequence_as_list[i][0] == '>':
                #get sequence length
                if i > 0:
                    sequences_dict[sequences_name] = one_line
                sequences_name = sequence_as_list[i]
                one_line = ""
                i += 2
            if i == (len(sequence_as_list) - 1):
                one_line += sequence_as_list[i].strip("\n")
                sequences_dict[sequences_name] = one_line

        files_dict[filename] = sequences_dict

    for motif in motifs_dict1:
        for filename in files_dict:
            motifs_dict1[motif][filename] = {}
            start_site_dict[filename] = {}
            for sequence in files_dict[filename]:
                motifs_dict1[motif][filename][sequence] = []
                start_site_dict[filename][sequence] = 0

    motifs_found1 = []
    highlighted_sequences = {}
    index_motif = 0
    keys = list(motif_possibilities.keys())
    for filename in files_dict:
        highlighted_sequences[filename] = {}
        for name in files_dict[filename]:
            sequence1 = files_dict[filename][name]
            highlighted_sequences[filename][name] = sequence1
            i = len(sequence1) - 1
            j = 0
            atgFound = False
            while i in range(len(sequence1)):
                if atgFound == True:
                    if (i - 4) >= 0:
                        potential_motif1 = sequence1[i - 4:i]
                        n_occurence = 0
                        for char in potential_motif1:
                            if "N" == char:
                                n_occurence += 1
                        if n_occurence < 2:
                            complement = GenerateComplement(potential_motif1)
                            for k in range(len(motif_possibilities.keys())):
                                if potential_motif1 in motif_possibilities[keys[k]] or complement in motif_possibilities[keys[k]]:
                                    motifs_dict1[keys[k]][filename][name].append(i)


                    i -= 1
                    j -= 1
                else:
                    potential_site = sequence1[(i - 2):(i + 1)]
                    if potential_site == "ATG" or potential_site == "CAT":
                        start_site = i - 2
                        start_site_dict[filename][name] = start_site
                        atgFound = True
                    if atgFound == False:
                        i -= 1
                    elif atgFound == True:
                        i = start_site

    cluster = {}
    window_in_list = {}
    for motif in motifs_dict1:
        cluster[motif] = {}
        window_in_list[motif] = {}
        print("For the motif " + motif + ":")
        for filename in files_dict:
            cluster[motif][filename] = {}
            window_in_list[motif][filename] = {}
            print("\tIn file " + filename + ":")
            for name in files_dict[filename]:
                print("\t\tIn sequence " + name.strip("\n") + ":")
                cluster[motif][filename][name] = {}
                window_in_list[motif][filename][name] = []
                sequence1 = files_dict[filename][name]
                motif_list = motifs_dict1[motif][filename][name]
                motif_index = 0
                cluster[motif][filename][name] = ""
                concatenated_windows = ""
                for index in motifs_dict1[motif][filename][name]:
                    downstream = start_site_dict[filename][name] - index
                    start = index
                    curr_window = sequence1[start:(index + window + 1)]
                    if motif_index < (len(motif_list) - 1):
                        upstream = (motif_list[motif_index] - motif_list[motif_index + 1]) - 4
                    else:
                        upstream = index - 4
                    not1st_downstream = (motif_list[motif_index - 1] - motif_list[motif_index]) - 4
                    #prevents first window from overlapping start site and also prevents overlap in upstream
                    if upstream >= (2*window) and downstream >= window and motif_index == 0: #done
                        curr_window = sequence1[start - window - 4:start - 4] + sequence1[start:start + window]
                        concatenated_windows += curr_window
                        window_in_list[motif][filename][name].append(curr_window)

                    elif upstream >= (2*window) and downstream < window and motif_index == 0:
                        #done
                        curr_window = sequence1[start - window - 4:start - 4] + sequence1[start:start + downstream]
                        concatenated_windows += curr_window
                        window_in_list[motif][filename][name].append(curr_window)

                    elif upstream < (2*window) and downstream >= window and motif_index == 0:
                        #done
                        if upstream > window:
                            beg_of_upstream_window = 0
                            if motif_index < (len(motif_list) - 1):
                                beg_of_upstream_window = start - 4 - (upstream - window)
                            else:
                                if (start - 4 - window) < 0:
                                    beg_of_upstream_window = 0
                                else:
                                    beg_of_upstream_window = start - window - 4
                            curr_window = sequence1[beg_of_upstream_window:start - 4] + sequence1[start:start + window]

                        else:
                            if motif_index < (len(motif_list) - 1):
                                beg_of_upstream_window = start - 4
                            else:
                                if (start - 4 - window) < 0:
                                    beg_of_upstream_window = 0
                                else:
                                    beg_of_upstream_window = start - window - 4
                            curr_window = sequence1[beg_of_upstream_window:start - 4] + sequence1[start:start + window]

                        concatenated_windows += curr_window
                        window_in_list[motif][filename][name].append(curr_window)

                    elif upstream < (2*window) and downstream < window and motif_index == 0:
                        #done
                        if upstream > window:
                            if motif_index < (len(motif_list) - 1):
                                beg_of_upstream_window = start - 4 - (upstream - window)

                            else:
                                if (start - 4 - window) < 0:
                                    beg_of_upstream_window = 0
                                else:
                                    beg_of_upstream_window = start - window - 4

                            curr_window = sequence1[beg_of_upstream_window:start - 4] + sequence1[start:start + downstream]
                        else:
                            if motif_index < (len(motif_list) - 1):
                                curr_window = sequence1[start:start + downstream]
                            else:
                                if (start - 4 - window) < 0:
                                    beg_of_upstream_window = 0
                                else:
                                    beg_of_upstream_window = start - window - 4
                                curr_window = sequence1[beg_of_upstream_window:start - 4] + sequence1[start:start + downstream]

                        concatenated_windows += curr_window
                        window_in_list[motif][filename][name].append(curr_window)

                    if motif_index > 0:
                        #if both windows are large enough to allow for normal window
                        if upstream >= (2*window) and not1st_downstream >= (2*window):
                            #done
                            curr_window = sequence1[start - window - 4:start - 4] + sequence1[start:start + window]
                            concatenated_windows += curr_window
                            #print(curr_window)
                            window_in_list[motif][filename][name].append(curr_window)
                            #CheckWindows(concatenated_windows, window, cluster)

                        #if window upstream is large enough to allow for a normal window to be taken but window downstream will lead to overlap unless adjusted
                        elif upstream >= (2*window) and not1st_downstream < (2*window):
                        #done
                            if not1st_downstream >= window:
                                curr_window = sequence1[start - window - 4:start - 4] + sequence1[start:start + window]

                            else:
                                if not1st_downstream > 0:
                                    curr_window = sequence1[start - window - 4:start - 4] + sequence1[start:start + not1st_downstream]
                                else:
                                    curr_window = sequence1[start - window - 4:start - 4]

                            concatenated_windows += curr_window
                            #print(curr_window)
                            window_in_list[motif][filename][name].append(curr_window)
                            #CheckWindows(concatenated_windows, window, cluster)

                        #if both windows up and downstream will lead to overlap unless adjusted
                        elif upstream < (2*window) and not1st_downstream < (2*window):
                            #done
                            if not1st_downstream >= window:
                                end_of_downstream_window = start + window

                            else:
                                if not1st_downstream > 0:
                                    end_of_downstream_window = start + not1st_downstream
                                else:
                                    end_of_downstream_window = start

                            if upstream > window:
                                #gets window upstream of size that allows for the next index to get downstream window of size of window variable without overlapping
                                if motif_index < (len(motif_list) - 1):
                                    beg_of_upstream_window = start - 4 - (upstream - window)
                                else:
                                    if (start - 4 - window) < 0:
                                        beg_of_upstream_window = 0
                                    else:
                                        beg_of_upstream_window = start - window

                            else:
                                if motif_index < (len(motif_list) - 1):
                                    beg_of_upstream_window = start - 4
                                else:
                                    if (start - 4 - window) < 0:
                                        beg_of_upstream_window = 0
                                    else:
                                        beg_of_upstream_window = start - window - 4

                            curr_window = sequence1[beg_of_upstream_window:start - 4] + sequence1[start:end_of_downstream_window]
                            concatenated_windows += curr_window
                            window_in_list[motif][filename][name].append(curr_window)

                        #if window downstream is large enough to allow for a normal window to be taken but window upstream will lead to overlap unless adjusted
                        elif upstream < (2*window) and not1st_downstream >= (2*window):
                            #done
                            if upstream > window:
                                if motif_index < (len(motif_list) - 1):
                                    beg_of_upstream_window = start - 4 - (upstream - window)
                                else:
                                    if (start - 4 - window) < 0:
                                        beg_of_upstream_window = 0
                                    else:
                                        beg_of_upstream_window = start - window - 4

                                curr_window = sequence1[beg_of_upstream_window:start - 4] + sequence1[start:start + window]

                            else:
                                if motif_index < (len(motif_list) - 1):
                                    curr_window = sequence1[start:start + window]
                                else:
                                    if (start - 4 - window) < 0:
                                        curr_window = sequence1[0:start - 4] + sequence1[start:start + window]
                                    else:
                                        curr_window = sequence1[start - window - 4:start - 4] + sequence1[start:start + window]
                            concatenated_windows += curr_window
                            window_in_list[motif][filename][name].append(curr_window)

                    motif_index += 1

                cluster[motif][filename][name] = concatenated_windows
                if len(window_in_list[motif][filename][name]) > 0:
                    co_motifs = FindRandMotifs(concatenated_windows, motif_size)
                    co_motifs1 = FindRandMotifs(window_in_list[motif][filename][name], motif_size)
                    print("List ")
                    for key in co_motifs1:
                        print(str(key) + ": " + str(co_motifs1[key]))
                    print("\nString ")
                    for key in co_motifs:
                        print(str(key) + ": " + str(co_motifs[key]))
'''
    for filename in files_dict:
        for name in files_dict[filename]:
            sequences = FindSeqWithMotifs("GATA", 2, "GGAA", 1, files_dict[filename][name], 30, start_site_dict[filename][name])
            print(sequences)
'''
#change function names to proper format "def word_word_word"

#gives index of co-motif occurences around a certain window of a motif

def ParseArgs(inputs):
    i = 0
    motifs_input = []
    onMotifs = False
    files_list = []
    onFiles = False
    window = 30
    onWindow = False
    motif_size = 4
    onSize = False
    while i < len(inputs):
        if inputs[i] == "-m":
            onMotifs = True
            onFiles = False
            onWindow = False
            onSize = False
            i += 1
        elif inputs[i] == "-f":
            onMotifs = False
            onFiles = True
            onWindow = False
            onSize = False
            i += 1
        elif inputs[i] == "-ws":
            onMotifs = False
            onFiles = False
            onWindow = True
            onSize = False
            i += 1
        elif inputs[i] == "-ms":
            onMotifs = False
            onFiles = False
            onWindow = False
            onSize = True
            i += 1

        if onMotifs == True:
            motifs_input.append(inputs[i])
        elif onFiles == True:
            files_list.append(inputs[i])
        elif onWindow == True:
            window = int(inputs[i])
        elif onSize == True:
            motif_size = int(inputs[i])

        i += 1
    return [motifs_input, files_list, window, motif_size]

#have not changed to be used as a independent function
def CheckWindows(concatenated_window, window, cluster):
    while k in range((index - window), (index + window + 1)):
        curr_potential_motif1 = sequence1[k:k + motif_size]
        complement = GenerateComplement(curr_potential_motif1)
        if curr_potential_motif1 in cluster[motif][filename][name][index]:
            cluster[motif][filename][name][index][curr_potential_motif1].append((index - k))
        elif complement in cluster[motif][filename][name][index]:
            cluster[motif][filename][name][index][complement].append((index - k))
        else:
            cluster[motif][filename][name][index][curr_potential_motif1] = [(index - k)]
        k += 1

#have not changed to be used as a independent function
def PrintHighlighted(motifs_dict1):

    for motif in motifs_dict1:
        for filename in files_dict:
            highlighted_sequences[filename] = {}
            for name in files_dict[filename]:
                print("\033[1;33m"+motif+"\033[3;37m")
                sequence = ""
                sequence = str(files_dict[filename][name])
                highlighted_motif = "\033[1;33m" + motif + "\033[3;37m"
                sequence = sequence.replace(motif, highlighted_motif)
                complement = GenerateComplement(motif)
                highlighted_complement = "\033[1;33m" + complement + "\033[3;37m"
                sequence = sequence.replace(complement, highlighted_complement)
                files_dict[filename][name] = sequence
                print(name)
                print(files_dict[filename][name])

    return sequence

def GenerateComplement(sequence):

    complement = ""
    i = len(sequence) - 1

    while i >= 0:
        if sequence[i] == "A":
            complement += "T"
        elif sequence[i] == "G":
            complement += "C"
        elif sequence[i] == "C":
            complement += "G"
        elif sequence[i] == "T":
            complement += "A"
        elif sequence[i] == "N":
            complement += "N"
        elif sequence[i] == "W":
            complement += "W"
        elif sequence[i] == "S":
            complement += "S"
        elif sequence[i] == "M":
            complement += "K"
        elif sequence[i] == "K":
            complement += "M"
        elif sequence[i] == "R":
            complement += "Y"
        elif sequence[i] == "Y":
            complement += "R"
        elif sequence[i] == "B":
            complement += "V"
        elif sequence[i] == "D":
            complement += "H"
        elif sequence[i] == "H":
            complement += "D"
        elif sequence[i] == "V":
            complement += "B"
        i -= 1

    return complement

def GenerateComplements(sequence):

    complements = [""]
    i = len(sequence) - 1

    while i < len(sequence) and i >= 0:
        h = 0
        if sequence[i] == "A":
            for j in range(len(complements)):
                complements[j] += "T"
        elif sequence[i] == "G":
            for k in range(len(complements)):
                complements[k] += "C"
        elif sequence[i] == "C":
            for m in range(len(complements)):
                complements[m] += "G"
        elif sequence[i] == "T":
            for p in range(len(complements)):
                complements[p] += "A"
        elif sequence[i] == "V":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements.insert((h), copy_seq)
                complements[h] += "G"
                complements[h + 1] += "T"
                complements[h + 2] += "C"
                h += 3
        elif sequence[i] == "H":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements.insert((h), copy_seq)
                complements[h] += "A"
                complements[h + 1] += "T"
                complements[h + 2] += "C"
                h += 3
        elif sequence[i] == "D":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements.insert((h), copy_seq)
                complements[h] += "A"
                complements[h + 1] += "T"
                complements[h + 2] += "G"
                h += 3
        elif sequence[i] == "B":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements.insert((h), copy_seq)
                complements[h] += "A"
                complements[h + 1] += "C"
                complements[h + 2] += "G"
                h += 3
        elif sequence[i] == "Y":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements[h] += "A"
                complements[h + 1] += "G"
                h += 2
        elif sequence[i] == "R":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements[h] += "C"
                complements[h + 1] += "T"
                h += 2
        elif sequence[i] == "K":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements[h] += "A"
                complements[h + 1] += "C"
                h += 2
        elif sequence[i] == "M":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements[h] += "G"
                complements[h + 1] += "T"
                h += 2
        elif sequence[i] == "S":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements[h] += "C"
                complements[h + 1] += "G"
                h += 2
        elif sequence[i] == "W":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements[h] += "A"
                complements[h + 1] += "T"
                h += 2
        elif sequence[i] == "N":
            while h < len(complements):
                copy_seq = complements[h]
                complements.insert((h), copy_seq)
                complements.insert((h), copy_seq)
                complements.insert((h), copy_seq)
                complements[h] += "A"
                complements[h + 1] += "T"
                complements[h + 2] += "C"
                complements[h + 3] += "G"
                h += 4
                #REMOVE FOR LOOPS
        i -= 1
    return complements

def FindRandMotifs(sequence, motif_size):

    most_common_motifs = {}
    motifs_found1 = {}
    motif_indexes = {}

    if type(sequence) is str:
        i = 0
        print(sequence)
        while i <= (len(sequence) - motif_size) and len(sequence) > 3:
            curr_potential_motif1 = sequence[i:i + motif_size]
            n_occurence = 0
            for char in curr_potential_motif1:
                if "N" == char:
                    n_occurence += 1
            if n_occurence < 2:
                complement = GenerateComplement(curr_potential_motif1)
                if curr_potential_motif1 in motifs_found1:
                    motifs_found1[curr_potential_motif1] += 1
                    motif_indexes[curr_potential_motif1].append(i)
                elif complement in motifs_found1:
                    motifs_found1[complement] += 1
                    motif_indexes[complement].append(i)
                else:
                    motifs_found1[curr_potential_motif1] = 1
                    motif_indexes[curr_potential_motif1] = [i]
            i += 1

        for j in range(10):
            key = max(motifs_found1.items(), key=operator.itemgetter(1))[0]
            motifs_found1.pop(key)
            most_common_motifs[key] = motif_indexes[key]

    if type(sequence) is list:
        index = 0
        for window in sequence:
            declared = False
            if len(window) > 3:
                sequence1 = window
                #print(window)
                i = 0
                #print(len(sequence1) - motif_size)
                while i <= (len(sequence1) - motif_size):
                    curr_potential_motif1 = sequence1[i:i + motif_size]
                    #print(curr_potential_motif1)
                    n_occurence = 0
                    for char in curr_potential_motif1:
                        if "N" == char:
                            n_occurence += 1
                    if n_occurence < 2:
                        complement = GenerateComplement(curr_potential_motif1)
                        if curr_potential_motif1 in motifs_found1:
                            motif_indexes[curr_potential_motif1].append(index)
                            motifs_found1[curr_potential_motif1] += 1
                        elif complement in motifs_found1:
                            motif_indexes[complement].append(index)
                            motifs_found1[complement] += 1
                        else:
                            motif_indexes[curr_potential_motif1] = [index]
                            motifs_found1[curr_potential_motif1] = 1

                    i += 1
                    index += 1
                index += motif_size - 1

        for j in range(10):
            key = max(motifs_found1.items(), key=operator.itemgetter(1))[0]
            motifs_found1.pop(key)
            most_common_motifs[key] = motif_indexes[key]

    return most_common_motifs


def FindSeqWithMotifs(motif1, num_motif1, motif2, num_motif2, sequence, window_to_search, start_site, motif_size):

    #list of sequences that fit the criteria
    list_sequences = []
    print(sequence)
    if start_site > 0:
        i = 0
        while i in range((start_site - window_to_search)):
            window = sequence[i:i + window_to_search]
            num_occurrences1 = 0
            num_occurrences2 = 0
            for j in range(len(window) - 3):
                potential_occurence = window[j:j + motif_size]
                if potential_occurence == motif1 or GenerateComplement(potential_occurence) == motif1:
                    num_occurrences1 += 1
                elif potential_occurence == motif2 or GenerateComplement(potential_occurence) == motif2:
                    num_occurrences2 += 1
            if num_occurrences1 == num_motif1 and num_occurrences2 == num_motif2:
                list_sequences.append(sequence[i:i + window_to_search])
            i += 1
    return list_sequences

main()

'''
>>CI-OTX 60

CGTTATCTCTAACGGAAGTTTTCGAAAAGGAAATTGTTCAATATCTAAGATAGGAATG

>>Co-Otx 60

CCGGTAATTACTCCATCAACAATAATTCAAGTTCATCTCAATTGAGGTAAGGAAAAACT
TCAATTTCTTAATTGAACAAATTTATTTAAAATCTAAATTGTTTATCAAGACGAGATGT
TTACTCTTTCGGTACACATTTATATTATAAACAAATTCCCCTTATGTTATGGGAACTAA
TATTCGAAGGCGTCTGGTCATTTATAGCAAGTAAAAAATACGTAAATAAATAAACAAAA
TACGAAGCCCGGAATAATAACATGCGATGAATGTATTCGGGTTTTGCTATTGTGTTACG
TAAAGAGCGCATTTTTTGTGTCAATAAATATTCTGCATTAGAATAAACGACAAAGCTAT
TGCACCATGGGAAAACTGCGTAGTAGGATTTCTTTTGGACAACGTAATAAGTTTGCGCA
ATTTTTCTTACAATTTGCAGCTCTAAAATATTTATAAAGACGGATATGAATTTAAAAAT
GATTTAAATTTATGAATTAATATATGGAATTAAGTTATTTCTAATAGATCATTTTCGGT
TGATAATATAATTCATATGAATTATTATATTAAATTTAGTTATAGTTAATAGATCATTT
TAGGATTTATAAATTAGTCTTAGCCATAAAAATATTCAGAAAATAATATTTGACCGCTG
GTCGCTCACGCATCACGTGATCACCATCCTGTTTAAATTATTGGACAAGATTCGTTTAA
GTTTACATAGGGTTCCTTGGTAGAGGATTAACTTTATTAATCCCTATAATTACCTTTGA
TTTCATTCGAATGTCGACCATTGTTCGGAAGCCAAGCAGAATGTTAACCAACGAAATCA
AAACCGTATTTTGACGCGTGTACTAATGCATAATTTGTCCATTGTGCTGCTTTATTTGC
TTTGGAGCCCAGAAAACGTGCACCGGTTTTTGCAGGCAGTAATTAGTCGCCGAAACAAA
AGTTCTAAAGTTCGTTTTCCGCACTCTTATAAACTACCTTGCAGCGCGCACATGTATGT
AATGAGGACAGCCGGTTATAACAATGCATCGGCAAGGTGCCGTATTTGAATAGGCCTTC
TCAGAGAGTGTGTGAAACTTTAAGTAATAGTTAATAAGAGTAAATCGTTAATATACGTA
TAAGCCGATCAGAAGGATATAACTTATTAAAATTATATAACTCTTTTTGTTTCAGGTTT
TAAGTTGAGTTTTTGTCAATTCACAAACGGAAATTATGTCTTATTTGAAGCCACCTCAT
TATGCCATG

>>Moccul.mesp.enhancer 60

AATCGGTGAAACCATATATTTCTATATATGTTGGGATTCGGTGGGATGAAATTATGTTA
TTTAGAACCGACTTACATTTGTGCATCTCTGACTGACTTGCACGACAAATATCGTTGTA
CCATCGGTTGGTTGATGGACGCCATCATCACCCAAATGATTGCACCACCGATCCCAAGT
GACCATGACGTGTGTTTGATGTATGGATCCGTGTCCCAACTGAAACACGAAATTGTTAA
AACGTAATATTAAACGTTGCAAATGGACAAGAGACTGAAAGTAGATGACAGCTTAAACT
AAATAATTTGTCAAGTGGCCGTTTGTTATAGTTTTGTATGTTTTGGTTAAATATTTCAC
TTCTTACTTGAAGAAATTGTCTCTTTTCCCAATTCTCATGGCCTCCAGAACTTGATCAA
AACCCCCAACAGTTGCACTTGTCTTTGCAATAACAAGAATCATTCCAAGGACCATTACA
CCAGCCTGAACAACATCTGTCCATATGACTCCTTTGATGCCGCCCTGTGAAGAAGAGGG
TTGAAACTATTGCATCAAACACATTAACAGAGAAATTCAACATTATACATCGTACCAAA
AAGTGTGTTGTGACAAAGATTTTAAAATTTTACAACTATAATACCGATCCGGAATGTTA
TTTCAGAAGAACCGAGGTACTGTCGGTAAGTGCTTGAAGGTCAAGAAATTGGCGAAGTA
TTGCACAGGTAGCAAAAGATGCAACATCTCGATTCTCGACGTAAAGTTGCAGACAATTA
ATTAATTAAACACTTCCGACATGATTGATCATCGCCGTTATCAGGAAAACTGATAACAT
TATCGCAAGTTTTATCGAGAACTGATTGCAAGGCGATAAGACTGGATTCCCTTCAACCA
AGATTGAAATAGCAGGAGGCAAAATGGGGTGAATTATGTCTGGACAAATGACAACTTTA
TGTACTATAAATCACGCAATTGTTTGTCATAAGCACACATAGTCAACGATAAAGTTATT
AATTCATTAAATTTAATCGTTTATATTATGAACAAACAAACGATTGTCATGAGAAGTTA
TTGCGAAGAAACAATCCAGTCTGTTCGATTGCAATGTCCAAATCAACAAACGAGCGTTC
CCGATGTATTTACTGACCGAGTTCCAAACAAACAAACCTTTATTGCGTATCCTGGGTCT
ACTGCTTCACATTTCCCTGC

>>Moocul.mesp.enhancer 60

NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTGTGTCGGATCAAGATCTCA
TCTCTACCCCAAACTCAAACCCCACCTACACTTGAGGTATTTGTCACAACGACATAACA
ACAAGAAATGGACGGATTGAGTGAAATTTCTAATGAATTTGTAATCAACTTATTCGTTT
TATTGACATCATTTGGCAACTGATGGTATATGGATCCATTATAGCACATCATACCATTA
CTTTATAATGGAATCCATACCATTAATGTTACTTGCGATGCGCACAGTGTAGCATATAC
ATAAGAATGGTTTCAGGCTATTCGACTATTGTACATACATTTGTGCGTCTCTCATTGAC
TTGCATGACAAATATCGTTGTACCATCGGTTGATTAATAGCTGCCATCATCACCCAGAT
GATGGCGCCACCGATCCCCAGTGACCAAGATGTGTGTTTGATGTATGGATCCGTATCCC
AACTGAAACATGAAATTGTGAAGCGTACTCAGTGCTGGGCAATAAACTGAAACGATGTG
AAAACATCCACTGCAAAATTTCAGTACAATTCGCATGGTTCTGATTAAGATAATGGAAT
ATTACAATTTCAGTAATAGAGTTGATAACGGAAAATAAGTATACTTATGGGTTTATACT
TATTTTATTATAAAGGTTGTTAAAAAAGTTAAAAAATAGTTGTGGATATCGGTCCAAAA
ACTTAAAATAACCCATCAAATAGCTTTCAATTCGTACTTGAAGAAGTTGTCTCTTTTTC
CAACCCTCATCGCCTCCATCACTTGATCAAAGCCTCCAACAATTGCACTTGTCTTCGCA
ATCACAAGAATCATTCCGAGTATCATGACGCCAGCCTGCACGACATCCGTCCATATGAC
TCCTTTGATGCCACCCTGTTGAAAGTGAAAGGTTATGTGATCGGTACGCAGATAAAGTA
TGGGAAATAAGACAGTCTTACATGACCGGTCACGGTCATGTAACCTAAACCATGACTTA
CACAGATCAAGATCTATTACAACATGCTATATACGGTTATAATATACATTGTGTTTCAA
TTGTCAATTAAGTACCACGAACCATCCATCACATTAAAGTTAATACTTTATAACAGAAG
CCAGCATTAATTCTGTCAAAAATATTTGCAATATTCACTGGTGGTTTCGATTTTAAGAT
ATCACAAATAATGTTTTCCAAAAATTTCTGTCAATGAAACCGAGATGTTGGATTACGTC
GGTAATTAGTAGCAGAAGTCGATGCAATATATCTCATCACATATGTTTACAGACAAATT
AATTAAACCCTTCCGACATGATTGGTCATCACCGCTATCACTAAAGCTGATAACATTAT
CGCAAGTTTTATCGAAAACTGATTGCAAGGCGATAAGACTGGATTGTCTTCGACGAAGC
TTGAAATAGCAGGAGGCAAAACGGAGTGAATTATGTCTGGACAAATACTGACTCTCTGT
TCTATAAATTCCTTCATTGTTTGTCATAAGCACACATTGTCAACGATAAAGTTATTAAT
TCATTAAATTTAATCGTTTATATTATGAACAAACAAACGATTGTCATGAGAAGTTATTA
TGAAGAAACAATCCAGTCAGTATCGTCACAATGTCCGAATCAACAAGCGGGCGTTCCCG
GTGTATTTACCACCCAAGTTCCAAACAAACAAACGTTTATTGCGTATC

>>Moocci.mesp.enhancer 60

GTGTTTAGCGATCACGGTTCAAGTTTAGGCTACATCGCTTTTTCCATCGGTTTGCCTAA
ACTTGCATGCTGTATACAGTATTAAGAAGGTCTTGTTTAACAATGTTCCTTAAATGTTA
TTGACATACAAAGCAAACTTTGTGTTAAAAGCCTACTTTAATGTGAAGACCAACATAAA
AAAATTATCGGATAGTATAACTTACTCAAAAAAGTTAAGTCTTTGTCCTCTTTGCAAAG
CTGCAGTAACTTTATCGAATCCCCCAAGATACACGGATGTTCTGATGATAACTAGAAGC
ATCCCCACTATCATTACTCCAGCTTGTACAACGTCTGTCCATATCACACCTTTTAAACC
ACCCTGGAAATTAAATCATTATTATTATTAAAGATAAAAACCATTTATAAGCCTATTTA
CAATACAGTTCCAAATTTAATCTGCCAATTCGATTGAACGTTTTTCAATGAAGTAAAAT
ATCTCTCACAATCTTGTATTGAGCCTTTAAAATAGCCTTTTAAATGAAAATCTATATCT
GACAAAATAGGCTTTGACGTGATTTGAACGCCGTACCTCATGCTCGTATAAATACCTAA
TTTGTTGCGGAATTTCTACTTCAATATAAAAAGTATGCAAATACATACATAATAGTTAT
TGAAAAGAATTTAAATCATCAAATTTTATGTCTGTTGTAGCCTAAATAGTGTTTTATTT
AAGTAAATCTACACTGCGCTTCAGTTCAAACCAACCGTGCCTTGAACTACGGGCAGTGC
GGTTGCTATAGACGTTATGCTTACGTAAACTAGTGCAAAGAATTTGTTGGAAAGTGTTA
ATTGATATCCGACCACCGACTGTGTTATCTTATCGTTTTGATAGATAAGAGCTGATAAC
ACTGATAATAGACTAATTAACTATAGCAGTTAATTTTTATTAAAATTTTCACTTTTCGG
ACATGAAAATTTGTCTAAAACCCAAACTTCCAAACTTCATCTATATAAGGAAAGAATTT
TATTGAAATGGCATACTCAATTCATCAGTATATATCATCTATATCTATATACGGTGGGT
AATACGTCGTAATTGATATCGGACAAAGTAAAAACTGCATATTGAAGCGTTAGAAGAAA
GTTGTTCAATTCAACATTTTAACAAAATG

>>Cirobu.mesp.enhancer 60

TTTTACATTTGAAATGTGATTAATTACGAAAATCCAGCGAATAGAATTGTCACAACAAG
TCATTAGCGACGGATATTTCGCCTTTGAAACTTAAAGGCGATAATGACTTTGCCCGTTT
CATGCGGCGATAAACGAACTAATTAGACACCTCCTACAGATATAATGGTAATTCAGAAT
CGTGTGGTTATGTAATTCACAAAAACATTTTAACAAAACAGGTTGATTTGAAACTTGTA
TTATGGAACTCAGTATGCATCAGGTAAGACGACGTGTACAATTTTTTCCAAAATAACGA
ATATA
'''
