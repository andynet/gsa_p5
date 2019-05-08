from bioparsers import read_fasta, read_fastq
import argparse
import json


def get_suffix_array(seq):

    suffix_pairs = sorted([(i, seq[i:]) for i in range(len(seq))], key=lambda suffix: suffix[1])

    suffix_array = [suffix_pairs[i][0] for i in range(len(seq))]
    suffices = [suffix_pairs[i][1] for i in range(len(seq))]

    return suffix_array, suffices


def get_c_table(seq):

    letter_count = dict()
    for i in range(len(seq)):
        if seq[i] not in letter_count.keys():
            letter_count[seq[i]] = 1
        else:
            letter_count[seq[i]] += 1

    sorted_letters = sorted(letter_count.items(), key=lambda letter: letter[0])
    c_table = dict()
    count = 0
    for i in range(len(sorted_letters)):
        c_table[sorted_letters[i][0]] = count
        count += sorted_letters[i][1]

        if sorted_letters[i][0] == '$':
            count -= 1

    return c_table


def get_o_table(seq, sa):

    bwt = ''.join([seq[i-1] for i in sa])

    o_table = dict()
    letters = set(seq)
    for letter in letters:
        o_table[letter] = []
        count = 0
        for i in range(len(bwt)):
            if bwt[i] == letter:
                count += 1
            o_table[letter].append(count)

    return o_table


def get_d_table(pattern, sa, c_table, o_prime_table):

    d_table = []

    L = 0
    R = len(sa) - 1
    ndiffs = 0

    for i in range(0, len(pattern)):

        try:
            letter = pattern[i]
            L = c_table[letter] + o_prime_table[letter][L-1] * (L != 0) + 1     # if L == 0: 0
            R = c_table[letter] + o_prime_table[letter][R]
        except KeyError:
            L = 0
            R = len(sa) - 1
            ndiffs += 1

        # if L > R: reset
        if L > R:
            L = 0
            R = len(sa) - 1
            ndiffs += 1

        d_table.append(ndiffs)

    return d_table


# def get_d_table_naive(pattern, string):
#
#     d_table = []
#
#     j = 0
#     ndiffs = 0
#
#     for i in range(len(pattern)):
#
#         if pattern[j:i+1] not in string:
#             j = i+1
#             ndiffs += 1
#
#         d_table.append(ndiffs)
#
#     return d_table


def inex_recur(pattern, sa, c_table, o_table, d_table, i, max_diffs, L, R, new_pattern, new_cigar):

    # print(f'{pattern} {i:3} {max_diffs:3} {L:3} {R:3} {new_pattern:10} {new_cigar:10}')
    intervals = []

    if i < 0:
        diffs = 0
    else:
        diffs = d_table[i]

    if max_diffs < diffs:  # if i = -1 then d_table[-1] == 0
        return intervals

    if i < 0:
        positions = sa[L:R+1]
        for position in positions:
            intervals.append((position, new_cigar))
        return intervals

    # deletion from pattern
    tmp_pattern = new_pattern[0:i] + new_pattern[i+1:]
    tmp_cigar = new_cigar[0:i] + 'D' + new_cigar[i+1:]
    result = inex_recur(pattern, sa, c_table, o_table, d_table, i-1, max_diffs-1, L, R, tmp_pattern, tmp_cigar)
    intervals = intervals + result

    letters = list(o_table.keys())[1:]
    for letter in letters:
        try:
            L_new = c_table[letter] + o_table[letter][L-1] * (L != 0) + 1
            R_new = c_table[letter] + o_table[letter][R]
        except KeyError:
            L_new = 1
            R_new = 0

        if L_new <= R_new:
            # insertion to pattern
            tmp_pattern = new_pattern[0:i] + letter + new_pattern[i:]
            tmp_cigar = new_cigar[0:i] + 'I' + new_cigar[i:]
            result = inex_recur(pattern, sa, c_table, o_table, d_table, i, max_diffs-1, L_new, R_new, tmp_pattern, tmp_cigar)
            intervals = intervals + result

            if letter == pattern[i]:
                # match
                tmp_pattern = new_pattern[0:i] + letter + new_pattern[i + 1:]
                tmp_cigar = new_cigar[0:i] + 'M' + new_cigar[i + 1:]
                result = inex_recur(pattern, sa, c_table, o_table, d_table, i-1, max_diffs, L_new, R_new, tmp_pattern, tmp_cigar)
                intervals = intervals + result
            else:
                # mismatch
                tmp_pattern = new_pattern[0:i] + letter + new_pattern[i + 1:]
                tmp_cigar = new_cigar[0:i] + 'M' + new_cigar[i + 1:]
                result = inex_recur(pattern, sa, c_table, o_table, d_table, i-1, max_diffs-1, L_new, R_new, tmp_pattern, tmp_cigar)
                intervals = intervals + result

    return intervals


def inexact_search(pattern, sa, c_table, o_table, o_prime_table, max_diffs):

    d_table = get_d_table(pattern, sa, c_table, o_prime_table)
    result = inex_recur(pattern, sa, c_table, o_table, d_table, len(pattern)-1, max_diffs, 0, len(sa)-1, pattern, 'M'*len(pattern))
    return result


def create_cigar(pseudo_cigar):

    cigar = ''

    last = pseudo_cigar[0]
    count = 1

    for i in range(1, len(pseudo_cigar)):

        if pseudo_cigar[i] != last:
            cigar += f'{count}{last}'
            last = pseudo_cigar[i]
            count = 1
        else:
            count += 1

    cigar += f'{count}{last}'

    return cigar


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('fasta', help='Fasta file')
    parser.add_argument('fastq', help='Fastq file')
    parser.add_argument('-p', help='preprocess', action='store_true')
    parser.add_argument('-d', help='distance', type=int)

    args = parser.parse_args()

    if args.p:
        fasta = read_fasta(args.fasta)

        sa_list, suff_list, c_list, o_list, o_prime_list = [], [], [], [], []
        for fasta_record in fasta:
            string = fasta_record.seq + '$'
            suffix_array, suffices = get_suffix_array(string)
            c_table = get_c_table(string)
            o_table = get_o_table(string, suffix_array)

            suffix_array_prime, suffices_prime = get_suffix_array(fasta_record.seq[::-1] + '$')
            o_prime_table = get_o_table(fasta_record.seq[::-1] + '$', suffix_array_prime)

            sa_list.append(suffix_array)
            suff_list.append(suffices)
            c_list.append(c_table)
            o_list.append(o_table)
            o_prime_list.append(o_prime_table)

        with open(args.fasta + '.sa', 'w') as outfile:
            json.dump(sa_list, outfile, sort_keys=True)
            print(f'Created file {outfile.name}')

        with open(args.fasta + '.c_tab', 'w') as outfile:
            json.dump(c_list, outfile, sort_keys=True)
            print(f'Created file {outfile.name}')

        with open(args.fasta + '.o_tab', 'w') as outfile:
            json.dump(o_list, outfile, sort_keys=True)
            print(f'Created file {outfile.name}')

        with open(args.fasta + '.o_prime_tab', 'w') as outfile:
            json.dump(o_prime_list, outfile, sort_keys=True)
            print(f'Created file {outfile.name}')

    else:
        with open(args.fasta + '.sa', 'r') as infile:
            sa_list = json.load(infile)

        with open(args.fasta + '.c_tab', 'r') as infile:
            c_list = json.load(infile)

        with open(args.fasta + '.o_tab', 'r') as infile:
            o_list = json.load(infile)

        with open(args.fasta + '.o_prime_tab', 'r') as infile:
            o_prime_list = json.load(infile)

        fasta = read_fasta(args.fasta)
        fastq = read_fastq(args.fastq)

        for i, fasta_record in enumerate(fasta):

            suffix_array = sa_list[i]
            c_table = c_list[i]
            o_table = o_list[i]
            o_prime_table = o_prime_list[i]

            for fastq_record in fastq:

                # print(f'{fastq_record.seq} in {fasta_record.seq} at:')

                results = inexact_search(fastq_record.seq, suffix_array, c_table, o_table, o_prime_table, args.d)
                # print(results)
                # result = set()
                # for position in list(positions):
                #     result = result.union(set(suffix_array[position[0]:position[1] + 1]))
                #
                # # print('result: ', result)
                for position, pseudo_cigar in results:
                    print(f'{fastq_record.sid}\t0\t{fasta_record.sid.strip()}\t{position + 1}\t0\t'
                          f'{create_cigar(pseudo_cigar)}\t*\t0\t0\t{fastq_record.seq}\t{fastq_record.qual}')


if __name__ == '__main__':
    main()