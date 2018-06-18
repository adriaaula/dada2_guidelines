import sys

def readFasta(fasta_file):
    """
    Input: fasta file
    Output: dict (header as key, seq as value)
    """
    fastaDict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            line = line.rstrip()
            if line.startswith('>'):
                header = line[1:] # remove '>'
                header = header.split(';')[0] # remove size annotations
                fastaDict[header] = ''
            else:
                fastaDict[header] += line
    return fastaDict

def readClusters(uc_file):
    """
    Read a UC format clusters file and return a dict
    with clustered seqID as key and seed seqID as value.
    Input: UC cluster file
    Output: clustersDict (key: clustered seqID, value: seed seqID)
    """
    clustersDict = {}
    with open(uc_file, 'r') as clusters:
        for line in clusters:
            line = line.strip()
            line = line.split('\t')
            clustered = line[8].split(';')[0] # remove size annotations
            seed = line[9].split(';')[0] # remove size annotations
            if line[0] == 'H':
                clustersDict[clustered] = seed
    return clustersDict

def main():
    """
    Takes fastaDict and clustersDict from previous functions
    and prints a data frame with ASV and OTU correspondence
    and their respective sequences
    """
    fasta_file = sys.argv[1]
    uc_file = sys.argv[2]
    fastaDict = readFasta(fasta_file)
    clustersDict = readClusters(uc_file)
    try:
        x = int(list(fastaDict.keys())[0])
        sortKey = int
    except ValueError:
        sortKey = str
    print('ASV', 'ASV.seq', 'OTU', 'OTU.seq', sep = '\t')
    for asv in sorted(fastaDict.keys(), key = sortKey):
        asvSeq = fastaDict[asv]
        if asv in clustersDict:
            otu = clustersDict[asv]
            otuSeq = fastaDict[otu]
        else:
            otu = asv
            otuSeq = asvSeq
        print(asv, asvSeq, otu, otuSeq, sep = '\t')

if __name__ == '__main__':
  main()