from pysqlite2 import dbapi2 as sqlite
import os
from astrolyze.setup.paths import prefix
# Keywords

connection = sqlite.connect(prefix.dataBase+'/parameter.db')
cursor = connection.cursor()

cursor.execute('CREATE TABLE Maps (id INTEGER PRIMARY KEY, Telescope VARCHAR(50), Species VARCHAR(50), uncertainty DOUBLE)')

telescope = 'IRAC'
species = '3.6MUM'
calError = 0.1
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))
species = '4.5MUM'
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))
species = '5.8MUM'
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))
species = '8MUM'
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))

telescope = 'MIPS'
species = '24MUM'
calError = 0.1
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))
species = '70MUM'
calError = 0.2
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))

telescope = 'PACS'
species = '100MUM'
calError = 0.2
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))
species = '160MUM'
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))

telescope = 'SPIRE'
species = '250MUM'
calError = 0.15
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))
species = '350MUM'
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))
species = '500MUM'
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))

telescope = 'GISMO'
species = '2MM'
calError = 0.15
cursor.execute('INSERT INTO Maps VALUES (null,?,?,?)',(telescope,species,calError))


connection.commit()

cursor.close()
connection.close()
