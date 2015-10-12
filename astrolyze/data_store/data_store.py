#!/usr/bin/env python

import fnmatch
import glob
import os

class ImportData(object):
    def __init__(self):
        pass

    def get_files_by_ending(self, directory, ending):
        """ Returns a list of all absolute path to the fits files"""

        file_list = [os.path.join(f)
                     for dirpath, dirnames, files in os.walk(directory)
                     for f in fnmatch.filter(files, '*.{}'.format(ending))]
        return file_list

if __name__ == "__main__":
    import_ = ImportData()
    directory = "/home/buchbend/Projects/Astronomy/"
    f = import_.get_files_by_ending(directory, "fits")
    print len(f)
    print len(set(f))
