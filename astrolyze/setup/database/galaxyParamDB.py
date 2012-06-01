from pysqlite2 import dbapi2 as sqlite
import os

from astrolyze.setup.paths import prefix

# Keywords

#Name
#MorphType
#Distance           [pc]
#VLSR               [km/s]
#Central Position    
#PA                 [deg]
#Inclination        [deg]
#R25                [min]

connection = sqlite.connect(prefix.dataBase+'/parameter.db')
cursor = connection.cursor()

cursor.execute('CREATE TABLE Galaxies (id INTEGER PRIMARY KEY,Name VARCHAR(50), MorphType VARCHAR(50), Distance DOUBLE, VLSR DOUBLE, Central Position VARCHAR(50), PA DOUBLE,Inclination FLOAT, R25 FLOAT)')

#M33 

cursor.execute('INSERT INTO Galaxies VALUES (null,?,?,?,?,?,?,?,?)',('M33','SA(s)cd',840e3,-179,'01:33:51.02 +30:39:36.7',-22.5,56,30.8))
#NGC3627

cursor.execute('INSERT INTO Galaxies VALUES (null,?,?,?,?,?,?,?,?)',('NGC3627','SAB(s)b',9.1e6,727,'11 20 15.027 +12 59 29.58',173,64,9.1))


connection.commit()

cursor.close()
connection.close()
