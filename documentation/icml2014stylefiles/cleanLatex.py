#! /usr/bin/python

# Python script to clean the temporary files from the LaTex build
# Author Mohanarajah Gajamohoan

import os

fileEndings = ['.out','*~','.synctex.gz','.bbl','.blg','.aux', '.log', '.dvi','.swp']
for root, dirs, files in os.walk('.'):
        for name in files:
                for endings in fileEndings:
                        if name.endswith(endings):
                                print('Removing ' + os.path.join(root, name) + ' ... '),
                                os.remove(os.path.join(root, name))
                                print('Done')

print 'Cleaning Done!'