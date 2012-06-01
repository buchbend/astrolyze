from pysqlite2 import dbapi2 as sqlite
import os
import constants as const
from astrolyze.setup.paths import prefix
# Keywords

# Name:  of the Line HAS TO BE UPPERCASE
# Frequency in GHz
# Wavelenght in m 


connection = sqlite.connect(prefix.dataBase+'/parameter.db')
cursor = connection.cursor()

cursor.execute('CREATE TABLE Lines (id INTEGER PRIMARY KEY,'
                                       'Name VARCHAR(50), '
                                        'Frequency FLOAT, '
                                        'Wavelenght Float)')

# Fill the table

hcop10 = 89.188523 * 1e9
cursor.execute('INSERT INTO Lines VALUES (null,?,?,?)',('HCOP10', hcop10, const.c/hcop10 ))

hcn10 = 88.6304156 * 1e9
cursor.execute('INSERT INTO Lines VALUES (null,?,?,?)',('HCN10', hcn10, const.c/hcn10))

CO1210 = 115.271204 * 1e9
cursor.execute('INSERT INTO Lines VALUES (null,?,?,?)',('12CO10', CO1210, const.c/CO1210))

CO1310 = 110.2013543 * 1e9
cursor.execute('INSERT INTO Lines VALUES (null,?,?,?)',('13CO10', CO1310, const.c/CO1310))

CO1221 = CO1210 * 2
cursor.execute('INSERT INTO Lines VALUES (null,?,?,?)',('12CO21', CO1221, const.c/CO1221))

HI =  21e-2 # 21cm line
cursor.execute('INSERT INTO Lines VALUES (null,?,?,?)',('HI', const.c/HI, HI))

Halpha =  457121.40*1e9 
cursor.execute('INSERT INTO Lines VALUES (null,?,?,?)',('HALPHA', Halpha, const.c/Halpha))


connection.commit()
cursor.close()
connection.close()
