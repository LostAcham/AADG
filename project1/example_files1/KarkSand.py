from array import array

def radixpass(char_indexes, sorted_indexes, key_values, char_indexes_length, alphabet_size):
    count_prefix = array("i", [0] * (alphabet_size + 1))
    for i in range(char_indexes_length):
        count_prefix[key_values[char_indexes[i]]] += 1

    prefix_sum = 0
    for i in range(alphabet_size + 1):
        freq, count_prefix[i] = count_prefix[i], prefix_sum
        prefix_sum += freq

    for i in range(char_indexes_length):
        sorted_indexes[count_prefix[key_values[char_indexes[i]]]] = char_indexes[i]
        count_prefix[key_values[char_indexes[i]]] += 1

def direct_kark_sort(word):
    alphabet = [None] + sorted(set(word))
    alphabet_size = len(alphabet)
    word_length = len(word)
    letter_dictionary = dict((char, index) for index, char in enumerate(alphabet))
    SA = array('i', [0] * (word_length + 3))
    kark_sort(array('i', [letter_dictionary[char] for char in word] + 3 * [0]), SA, word_length, alphabet_size)
    return SA[:word_length]

def kark_sort(char_indexes, SA, char_indexes_length, alphabet_size):
    nr_of_index_mod_0  = (char_indexes_length + 2) // 3
    nr_of_index_mod_1  = (char_indexes_length + 1) // 3
    nr_of_index_mod_2  = char_indexes_length // 3
    nr_of_index_mod_0_2 = nr_of_index_mod_0 + nr_of_index_mod_2
      
    SA_mod_1_2 = array('i', [0] * (nr_of_index_mod_0_2 + 3))
    SA_mod_0  = array('i', [0] * nr_of_index_mod_0)

    char_indexes_mod_1_2 = [i for i in range(char_indexes_length + (nr_of_index_mod_0 - nr_of_index_mod_1)) if i % 3] 
    char_indexes_mod_1_2.extend(3 * [0])
    char_indexes_mod_1_2 = array('i', char_indexes_mod_1_2)

    radixpass(char_indexes_mod_1_2, SA_mod_1_2, char_indexes[2:], nr_of_index_mod_0_2, alphabet_size)
    radixpass(SA_mod_1_2, char_indexes_mod_1_2, char_indexes[1:], nr_of_index_mod_0_2, alphabet_size)
    radixpass(char_indexes_mod_1_2, SA_mod_1_2, char_indexes, nr_of_index_mod_0_2, alphabet_size)

    nr_of_unique_trios = 0
    char_mod_0, char_mod_1, char_mod_2 = -1, -1, -1
    for suffix_index_mod_1_2 in range(nr_of_index_mod_0_2):
        if char_indexes[SA_mod_1_2[suffix_index_mod_1_2]] != char_mod_0 or char_indexes[SA_mod_1_2[suffix_index_mod_1_2] + 1] != char_mod_1 or char_indexes[SA_mod_1_2[suffix_index_mod_1_2] + 2] != char_mod_2 :
            nr_of_unique_trios += 1
            char_mod_0 = char_indexes[SA_mod_1_2[suffix_index_mod_1_2]]
            char_mod_1 = char_indexes[SA_mod_1_2[suffix_index_mod_1_2] + 1]
            char_mod_2 = char_indexes[SA_mod_1_2[suffix_index_mod_1_2] + 2]
        if SA_mod_1_2[suffix_index_mod_1_2] % 3 == 1:
            char_indexes_mod_1_2[SA_mod_1_2[suffix_index_mod_1_2] // 3] = nr_of_unique_trios
        else:
            char_indexes_mod_1_2[SA_mod_1_2[suffix_index_mod_1_2] // 3 + nr_of_index_mod_0] = nr_of_unique_trios

    if nr_of_unique_trios < nr_of_index_mod_0_2:
        kark_sort(char_indexes_mod_1_2, SA_mod_1_2, nr_of_index_mod_0_2, nr_of_unique_trios + 1)
        for suffix_index_mod_1_2 in range(nr_of_index_mod_0_2):
            char_indexes_mod_1_2[SA_mod_1_2[suffix_index_mod_1_2]] = suffix_index_mod_1_2 + 1
    else:
        for suffix_index_mod_1_2 in range(nr_of_index_mod_0_2):
            SA_mod_1_2[char_indexes_mod_1_2[suffix_index_mod_1_2] - 1] = suffix_index_mod_1_2

    char_indexes_mod_0 = array('i', [SA_mod_1_2[i] * 3 for i in range(nr_of_index_mod_0_2) if SA_mod_1_2[i] < nr_of_index_mod_0])
    radixpass(char_indexes_mod_0, SA_mod_0, char_indexes, nr_of_index_mod_0, alphabet_size)
  
    pointer_mod_0 = suffix_index_mod_0 = merged_index = 0
    pointer_mod_1_2 = nr_of_index_mod_0 - nr_of_index_mod_1
    while merged_index < char_indexes_length:
        suffix_index_mod_1_2 = SA_mod_1_2[pointer_mod_1_2] * 3 + 1 if SA_mod_1_2[pointer_mod_1_2] < nr_of_index_mod_0 else (SA_mod_1_2[pointer_mod_1_2] - nr_of_index_mod_0) * 3 + 2
        suffix_index_mod_0 = SA_mod_0[pointer_mod_0] if pointer_mod_0 < nr_of_index_mod_0 else 0

        if SA_mod_1_2[pointer_mod_1_2] < nr_of_index_mod_0 :
            suffix_mod_1_2_smaller = char_indexes_mod_1_2[SA_mod_1_2[pointer_mod_1_2] + nr_of_index_mod_0] <= char_indexes_mod_1_2[suffix_index_mod_0 // 3] if char_indexes[suffix_index_mod_1_2] == char_indexes[suffix_index_mod_0] else char_indexes[suffix_index_mod_1_2] < char_indexes[suffix_index_mod_0]
        elif char_indexes[suffix_index_mod_1_2] == char_indexes[suffix_index_mod_0]:
            suffix_mod_1_2_smaller = char_indexes_mod_1_2[SA_mod_1_2[pointer_mod_1_2] - nr_of_index_mod_0 + 1] <= char_indexes_mod_1_2[suffix_index_mod_0 // 3 + nr_of_index_mod_0] if char_indexes[suffix_index_mod_1_2 + 1] == char_indexes[suffix_index_mod_0 + 1] else char_indexes[suffix_index_mod_1_2 + 1] < char_indexes[suffix_index_mod_0 + 1]
        else:
            suffix_mod_1_2_smaller = char_indexes[suffix_index_mod_1_2] < char_indexes[suffix_index_mod_0]

    if suffix_mod_1_2_smaller:
        SA[merged_index] = suffix_index_mod_1_2
        pointer_mod_1_2 += 1
        if pointer_mod_1_2 == nr_of_index_mod_0_2:
            merged_index += 1
            while pointer_mod_0 < nr_of_index_mod_0:
                SA[merged_index] = SA_mod_0[pointer_mod_0]
                pointer_mod_0 += 1
                merged_index += 1
        
    else: 
        SA[merged_index] = suffix_index_mod_0
        pointer_mod_0 += 1
        if pointer_mod_0 == nr_of_index_mod_0:
            merged_index += 1
            while pointer_mod_1_2 < nr_of_index_mod_0_2:
                SA[merged_index] = (SA_mod_1_2[pointer_mod_1_2] * 3) + 1 if SA_mod_1_2[pointer_mod_1_2] < nr_of_index_mod_0 else ((SA_mod_1_2[pointer_mod_1_2] - nr_of_index_mod_0) * 3) + 2
                pointer_mod_1_2 += 1
                merged_index += 1
        merged_index += 1