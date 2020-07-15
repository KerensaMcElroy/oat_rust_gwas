#!/usr/bin/env python

# Copyright 2020
# Kerensa McElroy
# kerensa.mcelroy@csiro.au

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import re

def get_fasta(ref):
    fasta_dict={}
    with open(ref) as in_file:
        contig=''
        seq=''
        for i in in_file:
            if i[0]=='>':
                fasta_dict[contig]=seq
                contig=i
                seq=''
            else:
                seq+=i
        fasta_dict[contig]=seq
    return fasta_dict

def exclude(contigs, pattern):
    filtered=[]
    p = re.compile(pattern)
    for i in contigs:
        if p.search(i):
            pass
        else:
            filtered.append(i)

    return filtered

def fasta_out(filtered, fasta_dict, out):
    with open(out, 'w') as out_file:
        for f in filtered:
            seq=fasta_dict[f]
            out_file.write(f)
            out_file.write(seq)

def main(args):
    '''Main logic of program.'''
    
    in_fasta=get_fasta(args.infile)
    print(args.regex)
    filtered=exclude(in_fasta.keys(), args.regex)
    print(len(in_fasta.keys()), len(filtered))
    fasta_out(filtered, in_fasta, args.out) 

if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser(description = "Selects individual fasta records from a fasta \
                                                    file based on regular expression in record name.")

    parser.add_argument("-infile", help="fasta input file")
    parser.add_argument("-out", help="filtered output filename")
    parser.add_argument("-regex", help="python regular expression for filtering")

    args = parser.parse_args()

    main(args)
