#!/usr/bin/env python
import ConfigParser
import subprocess
import hashlib
import fnmatch
import glob
import os

from generaltools import file_tools
from generaltools import database_tools

USER = os.getenv("USER")
HOME = os.getenv("HOME")

class DataStore(object):
    def __init__(self):
        pass


class ImportData(database_tools.GenDb):
    def __init__(self):
        self.db = database_tools.SQLLiteConnection(
            "{}/.astrolyze/database/data_store.db".format(HOME)
        )

    def setup_map_store_table(self):
        """Intialize table for proceess coordination

        The central table that is used to store information about the maps has
        the following entries:

            - path
            - map name
            - hash
            - parent
        """
        table_name = "map_store"
        sql = "CREATE TABLE if not exists {0}"\
              "(id INTEGER PRIMARY KEY NOT NULL AUTO_INCREMENT, "\
              "full_path STRING, "\
              "map_name TEXT, "\
              "hash INTEGER, "
              "parent INT," \
              ")"\
              "".format(table_name, self.sources_table,
                        self.obsmodes_table)
        self.execute(sql)
        self.commit()



if __name__ == "__main__":
    import_ = ImportData()
    directory = "/home/buchbend/Projects/Astronomy/"
    files_ = file_tools.find_files_by_ending_in_directory(directory=directory,
                                                          ending="fits")
    for file_ in set(files_):
        subprocess.call("cp {} ~/data_store".format(file_), shell=True)
