#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 18:05:31 2019

@author: julienbryios
"""

import os
import sys
import hashlib
import send2trash

def hashfile(path, blocksize = 65536):
    afile = open(path, 'rb')
    hasher = hashlib.md5()
    buf = afile.read(blocksize)
    while len(buf) > 0:
        hasher.update(buf)
        buf = afile.read(blocksize)
    afile.close()
    return hasher.hexdigest()

def findDup(parentFolder):
    # Dups in format {hash:[names]}
    dups = {}
    for dirName, subdirs, fileList in os.walk(parentFolder):
        print('Scanning %s...' % dirName)
        for filename in fileList:
            # Get the path to the file
            path = os.path.join(dirName, filename)
            # Calculate hash
            file_hash = hashfile(path)
            # Add or append the file path
            if file_hash in dups:
                print(path)
                send2trash.send2trash(path)
            else:
                dups[file_hash] = [path]
    return dups

if __name__ == "__main__":
    findDup('/Users/julienbryios/Dropbox/Gabrielle/')