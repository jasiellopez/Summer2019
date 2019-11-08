import sys
import operator
import collections

def main():
    #the goal of this script is to identify the most commonly-occurring DNA sequences (motifs) 
    #within a certain range (given in terms of base pairs) of a known co-motif (also given)
    inputs = ParseArgs(sys.argv)
    motifs_input = inputs[0]
    files_list = inputs[1]
    window = inputs[2]
    motif_size = inputs[3]

    motif_dict = {}
    motif_possibilities = {}
    start_site_dict = {}
    
    for motif in motifs_input:
        motif_possibilities[motif] = GenerateComplements(motif)
        motif_dict[motif] = {}

    sequences_dict = {} #dictionary where key is sequence name and value is the sequence
    files_dict = {} #dictionary where the key is the file and value is a dictionary of its sequences

    for filename in files_list:
        parseSeqs(filename, sequences_dict)
        files_dict[filename] = sequences_dict
        sequences_dict = {} #reset sequence dictionary for a new file
        
    #initialize dictionaries which contain 
    for motif in motif_dict:
        for filename in files_dict:
            motif_dict[motif][filename] = {}
            start_site_dict[filename] = {}
            for sequence in files_dict[filename]:
                #motif_dict will contain index where motif occurs
                motif_dict[motif][filename][sequence] = []
                #start_site_dict contains the predicted start site for every sequence in every file
                start_site_dict[filename][sequence] = 0
                
    #find where each motif occurs in each sequence
    for filename in files_dict:
        for name in files_dict[filename]:
            findMotifsInSeqs(motif_dict, motif_possibilities, start_site_dict, files_dict, filename, name)
            
    #this entire loop creates a concatenated sequence which created from the 
    #sequences n base pairs up and downstream of a known "motif" occurrence.
    #overlap is accounted for and if this occurs, the algorithm pulls less 
    #than "window" base pairs
    window_as_string = {}
    window_in_list = {}
    for motif in motif_dict:
        window_as_string[motif] = {}
        window_in_list[motif] = {}
        print("For the motif " + motif + ":")
        for filename in files_dict:
            window_as_string[motif][filename] = {}
            window_in_list[motif][filename] = {}
            print("\tIn file " + filename + ":")
            for name in files_dict[filename]:
                print("\t\tIn sequence " + name.strip("\n") + ":")
                sequence1 = files_dict[filename][name]
                motif_list = motif_dict[motif][filename][name]
                #grab windows up and down stream of motif and concatenates them
                createSubStrings(start_site_dict, window_as_string, window_in_list, sequence1, motif_list, filename, name, motif, window)
                if len(window_in_list[motif][filename][name]) > 0:
                    #finds 10 most common motifs of a previously specified size, one using a string 
                    #and another using a list (in order to see how many accidental motifs were found in string case)
                    co_motifs = FindRandMotifs(window_as_string[motif][filename][name], motif_size)
                    co_motifs1 = FindRandMotifs(window_in_list[motif][filename][name], motif_size)
                    print("List ")
                    for key in co_motifs1:
                        print(str(key) + ": " + str(co_motifs1[key]))
                    print("\nString ")
                    for key in co_motifs:
                        print(str(key) + ": " + str(co_motifs[key]))


def createSubStrings(start_site_dict, window_as_string, window_in_list, sequence1, motif_list, filename, name, motif, window):

    motif_index = 0
    window_in_list[motif][filename][name] = []
    window_as_string[motif][filename][name] = ""
    concatenated_windows = ""
    for index in motif_list:
        curr_motif_size = len(motif)
        downstream = start_site_dict[filename][name] - index
        start = index
        curr_window = sequence1[start:(index + window + 1)]
        curr_window1 = ""
        if motif_index < (len(motif_list) - 1):
            upstream = (motif_list[motif_index] - motif_list[motif_index + 1]) - curr_motif_size
        else:
            upstream = index - curr_motif_size
        not1st_downstream = (motif_list[motif_index - 1] - motif_list[motif_index]) - curr_motif_size
        #prevents first window from overlapping start site and also prevents overlap in upstream
        if upstream >= (2*window) and downstream >= window and motif_index == 0:
            curr_window = sequence1[start - window - curr_motif_size:start - curr_motif_size]
            curr_window1 = sequence1[start:start + window]
            concatenated_windows += curr_window + curr_window1
            window_in_list[motif][filename][name].append(curr_window)
            window_in_list[motif][filename][name].append(curr_window1)
        #accounts for smaller than "window" sequence size downstream
        elif upstream >= (2*window) and downstream < window and motif_index == 0:
            curr_window = sequence1[start - window - curr_motif_size:start - curr_motif_size]
            curr_window1 = sequence1[start:start + downstream]
            concatenated_windows += curr_window + curr_window1
            window_in_list[motif][filename][name].append(curr_window)
            window_in_list[motif][filename][name].append(curr_window1)
        #accounts for smaller than "window" sequence size upstream
        elif upstream < (2*window) and downstream >= window and motif_index == 0:
            if upstream > window:
                beg_of_upstream_window = 0
                if motif_index < (len(motif_list) - 1):
                    beg_of_upstream_window = start - curr_motif_size - (upstream - window)
                else:
                    if (start - curr_motif_size - window) < 0:
                        beg_of_upstream_window = 0
                    else:
                        beg_of_upstream_window = start - window - curr_motif_size
                curr_window = sequence1[beg_of_upstream_window:start - curr_motif_size]
                curr_window1 = sequence1[start:start + window]

            else:
                if motif_index < (len(motif_list) - 1):
                    beg_of_upstream_window = start - curr_motif_size
                else:
                    if (start - curr_motif_size - window) < 0:
                        beg_of_upstream_window = 0
                    else:
                        beg_of_upstream_window = start - window - curr_motif_size
                curr_window = sequence1[beg_of_upstream_window:start - curr_motif_size]
                curr_window1 = sequence1[start:start + window]

            concatenated_windows += curr_window + curr_window1
            window_in_list[motif][filename][name].append(curr_window)
            window_in_list[motif][filename][name].append(curr_window1)
        #accounts for case when both up and downstream sequences are less than n base pairs
        elif upstream < (2*window) and downstream < window and motif_index == 0:
            if upstream > window:
                if motif_index < (len(motif_list) - 1):
                    beg_of_upstream_window = start - curr_motif_size - (upstream - window)

                else:
                    if (start - curr_motif_size - window) < 0:
                        beg_of_upstream_window = 0
                    else:
                        beg_of_upstream_window = start - window - curr_motif_size

                curr_window = sequence1[beg_of_upstream_window:start - curr_motif_size]
                curr_window1 = sequence1[start:start + downstream]
                concatenated_windows += curr_window + curr_window1
                window_in_list[motif][filename][name].append(curr_window)
                window_in_list[motif][filename][name].append(curr_window1)
            else:
                if motif_index < (len(motif_list) - 1):
                    curr_window = sequence1[start:start + downstream]
                else:
                    if (start - curr_motif_size - window) < 0:
                        beg_of_upstream_window = 0
                    else:
                        beg_of_upstream_window = start - window - curr_motif_size
                    curr_window = sequence1[beg_of_upstream_window:start - curr_motif_size]
                    curr_window1 = sequence1[start:start + downstream]

            concatenated_windows += curr_window + curr_window1
            window_in_list[motif][filename][name].append(curr_window)
            window_in_list[motif][filename][name].append(curr_window1)

        if motif_index > 0:
            #if both windows are large enough to allow for normal window
            if upstream >= (2*window) and not1st_downstream >= (2*window):
                curr_window = sequence1[start - window - curr_motif_size:start - curr_motif_size]
                curr_window1 = sequence1[start:start + window]
                concatenated_windows += curr_window + curr_window1
                window_in_list[motif][filename][name].append(curr_window)
                window_in_list[motif][filename][name].append(curr_window1)

            #if window upstream is large enough to allow for a normal window to be taken but window downstream will lead to overlap unless adjusted
            elif upstream >= (2*window) and not1st_downstream < (2*window):
                if not1st_downstream >= window:
                    curr_window = sequence1[start - window - curr_motif_size:start - curr_motif_size] + sequence1[start:start + window]

                else:
                    if not1st_downstream > 0:
                        curr_window = sequence1[start - window - curr_motif_size:start - curr_motif_size]
                        curr_window1 = sequence1[start:start + not1st_downstream]
                    else:
                        curr_window = sequence1[start - window - curr_motif_size:start - curr_motif_size]
                        curr_window1 = ""
                concatenated_windows += curr_window + curr_window1
                window_in_list[motif][filename][name].append(curr_window)
                window_in_list[motif][filename][name].append(curr_window1)

            #if both windows up and downstream will lead to overlap unless adjusted
            elif upstream < (2*window) and not1st_downstream < (2*window):
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
                        beg_of_upstream_window = start - curr_motif_size - (upstream - window)
                    else:
                        if (start - curr_motif_size - window) < 0:
                            beg_of_upstream_window = 0
                        else:
                            beg_of_upstream_window = start - window

                else:
                    if motif_index < (len(motif_list) - 1):
                        beg_of_upstream_window = start - curr_motif_size
                    else:
                        if (start - curr_motif_size - window) < 0:
                            beg_of_upstream_window = 0
                        else:
                            beg_of_upstream_window = start - window - curr_motif_size

                curr_window = sequence1[beg_of_upstream_window:start - curr_motif_size]
                curr_window1 = sequence1[start:end_of_downstream_window]
                concatenated_windows += curr_window + curr_window1
                window_in_list[motif][filename][name].append(curr_window)
                window_in_list[motif][filename][name].append(curr_window1)

            #if window downstream is large enough to allow for a normal window to be taken but window upstream will lead to overlap unless adjusted
            elif upstream < (2*window) and not1st_downstream >= (2*window):
                if upstream > window:
                    if motif_index < (len(motif_list) - 1):
                        beg_of_upstream_window = start - curr_motif_size - (upstream - window)
                    else:
                        if (start - curr_motif_size - window) < 0:
                            beg_of_upstream_window = 0
                        else:
                            beg_of_upstream_window = start - window - curr_motif_size

                    curr_window = sequence1[beg_of_upstream_window:start - curr_motif_size]
                    curr_window1 = sequence1[start:start + window]

                else:
                    if motif_index < (len(motif_list) - 1):
                        curr_window = sequence1[start:start + window]
                        curr_window1 = ""
                    else:
                        if (start - curr_motif_size - window) < 0:
                            curr_window = sequence1[0:start - curr_motif_size]
                            curr_window1 = sequence1[start:start + window]
                        else:
                            curr_window = sequence1[start - window - curr_motif_size:start - curr_motif_size]
                            curr_window1 = sequence1[start:start + window]
                concatenated_windows += curr_window + curr_window1
                window_in_list[motif][filename][name].append(curr_window)
                window_in_list[motif][filename][name].append(curr_window1)

        motif_index += 1

    window_as_string[motif][filename][name] = concatenated_windows
    

def findMotifsInSeqs(motif_dict, motif_possibilities, start_site_dict, files_dict, filename, name):
    atgFound = False
    for motif in motif_dict:
        sequence1 = files_dict[filename][name]
        i = len(sequence1) - 1
        j = 0
        while i in range(len(sequence1)):
            if atgFound == True:
                if (i - len(motif)) >= 0:
                    potential_motif1 = sequence1[i - len(motif):i]
                    complement = GenerateComplement(potential_motif1)
                    if potential_motif1 in motif_possibilities[motif] or complement in motif_possibilities[motif]:
                        motif_dict[motif][filename][name].append(i)

                i -= 1
                j -= 1

            else:
                #find start site (assuming first ATG is start site and that DNA is 5'-3') if hasn't been found
                potential_site = sequence1[(i - 2):(i + 1)]
                if potential_site == "ATG" or potential_site == "CAT":
                    start_site = i - 2
                    start_site_dict[filename][name] = start_site
                    atgFound = True
                if atgFound == False:
                    i -= 1
                elif atgFound == True:
                    i = start_site

def parseSeqs(filename, sequences_dict):
    file = open(filename, "r")
    sequence_as_list = file.readlines()
    one_line = ""
    i = 0
    #parses sequences in every file. Adds the sequence name as the key in dictionary 
    #"sequences_dict" and value is the sequence itself
    while i < len(sequence_as_list):
        if len(sequence_as_list[i]) > 0 and sequence_as_list[i][0] != '>':
            one_line += sequence_as_list[i].strip("\n")
            i += 1
        elif sequence_as_list[i][0] == '>':
            if i > 0:
                sequences_dict[sequences_name] = one_line
            sequences_name = sequence_as_list[i]
            one_line = ""
            i += 2
        if i == (len(sequence_as_list) - 1):
            one_line += sequence_as_list[i].strip("\n")
            sequences_dict[sequences_name] = one_line

#parses arguments given at program call
def ParseArgs(inputs):
    i = 0
    if len(inputs):
        print("Hello, there are 4 parameters available to you: known motifs (should be preceded by -m), files (please input full path and precede by -f), window size to search around motifs (-ws), and size which you want co-motifs to be (-ms).")

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
        elif inputs[i] == "-help":
            print("Hello, there are 4 parameters available to you: known motifs (should be preceded by -m), files (please input full path and precede by -f), window size to search around motifs (-ws), and size which you want co-motifs to be (-ms).")

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

#have not changed to be used as an independent function
"""
def CheckWindows(concatenated_window, window, cluster, index):
    k = 0
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
"""

"""
#have not changed to be used as an independent function, 
#it is an algorithm to print out found sequence with highlighted motifs to terminal
def PrintHighlighted(motif_dict):

    for motif in motif_dict:
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
"""

#create complement strand of sequence given
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

#converts the non-ATGC characters to ATGC based on universal DNA code
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
        i -= 1
    return complements

#finds motif candidates of "motif_size" in sequence given to function and 
#return 10 most common motifs found
def FindRandMotifs(sequence, motif_size):

    most_common_motifs = {}
    motifs_found1 = {}
    motif_indexes = {}

    if type(sequence) is str:
        i = 0
        j = 0
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
            if len(window) > 3:
                sequence1 = window
                i = 0
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

#finds sequences of DNA of size "window_to_search" which contain a "num_motif1" amount of 
#"motif1" occurrences and a "num_motif2" amount of "motif2" occurences
def FindSeqWithMotifs(motif1, num_motif1, motif2, num_motif2, sequence, window_to_search, start_site, motif_size, motif_possibilities):

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
                if potential_occurence in motif_possibilities[motif1] or GenerateComplement(potential_occurence) in motif_possibilities[motif1]:
                    num_occurrences1 += 1
                elif potential_occurence in motif_possibilities[motif2] or GenerateComplement(potential_occurence) in motif_possibilities[motif2]:
                    num_occurrences2 += 1
            if num_occurrences1 == num_motif1 and num_occurrences2 == num_motif2:
                list_sequences.append(sequence[i:i + window_to_search])
            i += 1
    return list_sequences

main()
