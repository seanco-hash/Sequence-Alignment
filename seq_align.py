import argparse
from itertools import groupby
import pandas as pd
import AlignFactory


def fastaread(fasta_name):
    """
    Read all sequences from a given fasta file path
    :param fasta_name: The path to the fasta file
    :return: An iterator to the (header, seq) of the sequences present in the file
    """
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


def main():
    """
    This program made in order to find the optimal two sequences alignment.
    It gets two sequences in fasta files, score matrix and wanted alignment
    type and prints the best alignment of the wanted type and its score.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')

    # Holds the given arguments
    args = parser.parse_args()

    # extract the sequences and the scores table
    s = fastaread(args.seq_a).__next__()[1]
    t = fastaread(args.seq_b).__next__()[1]
    scores = pd.read_csv(args.score, delimiter='\t')

    # Get the wanted alignment type object from the factory
    alignm = AlignFactory.getAlign(s, t, args.align_type, scores)

    # print the best alignment and score
    alignm.doAlign()

if __name__ == '__main__':
    main()
