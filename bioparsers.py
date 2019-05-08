class FastaRecord:

    def __init__(self, sid, seq):
        self.sid = sid
        self.seq = seq


def read_fasta(file):

    fasta = []

    with open(file) as f:
        lines = f.readlines()

    sid = None
    seq = ''

    for line in lines:

        if line.startswith('>'):
            if sid is not None:
                fasta.append(FastaRecord(sid, seq))

            sid = line.strip().strip('>')
            seq = ''
        elif line.startswith(';'):
            continue
        else:
            seq += line.strip()

    fasta.append(FastaRecord(sid, seq))

    return fasta


class FastqRecord:

    def __init__(self, sid, seq, qual):
        self.sid = sid
        self.seq = seq
        self.qual = qual


def read_fastq(file):

    fastq = []

    with open(file) as f:
        lines = f.readlines()

    for i in range(0, len(lines), 4):

        sid = lines[i].strip().strip('@')
        seq = lines[i+1].strip()
        qual = lines[i+3].strip()

        fastq.append(FastqRecord(sid, seq, qual))

    return fastq
