"""
Parse text-formatted(nofmt) blast output files to generate blast hit
sequence as csv files

"""

import sys
import glob

from Bio.Blast import NCBIStandalone


def getNposition(query, sbjct):
    aster = []
    for t, pair in enumerate(zip(query, sbjct), 1):
        q, s = pair
        if q == 'N':
            aster.append([t, q, s])
    return aster

if __name__ == "__main__":

    blast_parser = NCBIStandalone.BlastParser()

    GLB = glob.glob('*.nofmt')
    for glb in GLB:
        handle = open(glb, 'r')
        blast_iterator = NCBIStandalone.Iterator(handle, blast_parser)
        total = 0
        E_VALUE_THRESH = 1
        with open(glb[:-6] + ".csv", 'w') as outcsv:
            for blast_record in blast_iterator:
                print ("query:", blast_record.query)
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        if hsp.expect < E_VALUE_THRESH:
                            if hsp.align_length < 100:  
                                continue
                            if len(str(hsp.sbjct)) == 0:
                                continue
                            print(hsp.query)
                            print(hsp.match)
                            print(hsp.sbjct)
                            print('')

                            aster = getNposition(hsp.query, hsp.sbjct)
                            if len(aster) > 1:  # marker has SNP at multiple positions
                                for mult, q in enumerate(aster, 1):
                                    nuc = q[-1]
                                    outcsv.write(
                                        blast_record.query + '|' + str(mult) + ',' + nuc + '\n')
                                    mult += 1
                                    total += 1

                            else:
                                nuct = [q[-1] for q in aster][0]
                                outcsv.write(
                                    blast_record.query + ',' + nuct + '\n')
                                total += 1
        print ("Done Written:", total)
