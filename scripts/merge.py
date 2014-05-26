#!/usr/bin/env python3

import sys
import os.path as path
import itertools

REMAINING_KEYS = ['STANDARD', 'COLUMNID', 'CROSS', 'CROSSED']
REMOVED_KEYS = ['IGNORE', 'COMPLEX']

def usage():
    print("""\
Usage:
  python3 merge.py <filename>...

  Simplify PAD files by applying the information in the
  IGNORE and COMPLEX meta rows to the data rows.
""")


def remove_columns(data, specials):
    toIgnore = []
    for i, cell in enumerate(specials['IGNORE']):
        if not cell.strip():
            print("empty cells in IGNORE")
        if cell.strip() and cell in 'xX':
            toIgnore.append(i)
    toIgnore.reverse()
    extra = []
    for key in specials:
        if key != 'COLUMNID':
            extra.append(specials[key])
    for l in itertools.chain(data, extra):
        for i in toIgnore:
            try:
                del l[i]
            except IndexError:
                print('Row too short: %i %s' % (i, l))
    return len(toIgnore)


def merge_complex_cols(data, specials):
    ranges = []
    for i, cell in enumerate(specials['COMPLEX']):
        if cell == '<':
            ranges.append(i)
        elif cell == '>':
            try:
                start = ranges.pop()
                if not isinstance(start, int):
                    raise IndexError
                ranges.append((start, i+1))
            except IndexError:
                print("Inconsistent COMPLEX row in file '%s'" % name)
                sys.exit(1)
    ranges.reverse()
    #special rows
    for key in specials:
        if key in REMOVED_KEYS or key == 'COLUMNID':
            continue
        row = specials[key]
        for x, y in ranges:
            replacement = ''.join(row[x:y]).replace('.', '').replace('-', '')
            #special problem with CROSS(ED) meta rows: avoid '++' entries
            if replacement.count('+') == len(replacement) and \
              len(replacement) > 1:
              replacement = '+'
            row[x:y] = [replacement]
    #other rows
    for row in data:
        for x, y in ranges:
            snip = row[x:y]
            for i in range(len(snip)-1, -1, -1):
                if snip[i] == '-':
                    del snip[i]
            if len(snip) == 0:
                snip.append('-')
            row[x:y] = [''.join(snip)]
    
    return sum(y-x-1 for x, y in ranges)
            

def merge(name):
    print("processing " + name)
    sourceName = path.abspath(name)
    destName = path.abspath(path.basename(name))
    if sourceName == destName:
        print("Error: Will not overwrite source file '%s'" % name)
        sys.exit(1)
    #read in
    source = open(sourceName, 'r', encoding='utf8')
    lines = []
    data = []
    specials = {}
    for l in source.readlines():
        if l.startswith('#'):
            lines.append(l)
            continue
        l = [x.strip() for x in l.split('\t')]
        if l[0] == ':ANN':
            key = l[1].rstrip('.')
            specials[key] = l
            if key in REMAINING_KEYS:
                lines.append(l)
        else:
            data.append(l)
            lines.append(l)
    source.close()
    #modify
    reduction = 0
    if 'IGNORE' in specials:
        reduction += remove_columns(data, specials)
    if 'COMPLEX' in specials:
        reduction += merge_complex_cols(data, specials)
    if reduction:
        specials['COLUMNID'][-reduction:] = []
    #write out
    out = open(destName, 'w', encoding='utf8')
    for l in lines:
        if isinstance(l, list):
            out.write('\t'.join(l) + '\n')            
        else:
            out.write(l)
    out.close()


def main(args):
    if not args or args[0] in ['-h', '--help']:
        usage()
        sys.exit(1)
    for name in args:
        merge(name)


if __name__ == '__main__':
    main(sys.argv[1:])
