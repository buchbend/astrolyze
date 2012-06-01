from pysqlite2 import dbapi2 as sqlite
from astrolyze.setup.paths import init_prefix

print sysPrefix

class galaxyParams:
    def __init__(self, name):
        connection = sqlite.connect(str(init_prefix.sysPrefix)+'/HERMES/dataBase/galaxyParameter.db')
        cursor = connection.cursor()
        cursor.execute("SELECT * FROM Galaxies WHERE Name = ?",(name,))
        self.params =  cursor.fetchall()[0]
        self.name = self.params[1]
        self.type = self.params[2]
        self.distance = self.params[3]
        self.vlsr = self.params[4]
        self.centralPosition = self.params[5]
        self.pa = self.params[6]
        self.inclination = self.params[7]
        self.R25 =self.params[8]
        

    def showParams(self):
        
        print 'Name = '+str(self.name)
        print 'Type = '+str(self.type)
        print 'Distance = '+str(self.distance)
        print 'VLSR = '+str(self.vlsr)
        print 'Central Position = '+str(self.centralPosition)
        print 'PA = '+str(self.pa)
        print 'Inclination = '+str(self.inclination)
        print 'R25 = '+str(self.R25)


class lineParameter:
    def __init__(self, name):
        connection = sqlite.connect(str(init_prefix.sysPrefix)+'/HERMES/dataBase/lineParameter.db')
        cursor = connection.cursor()
        cursor.execute("SELECT * FROM Lines WHERE Name = ?",(name,))
        self.params =  cursor.fetchall()[0]
        self.name = self.params[1]
        self.frequency = self.params[2]
        self.wavelenght = self.params[3]

class lineParameter:
    def __init__(self, name):
        connection = sqlite.connect(str(init_prefix.sysPrefix)+'/HERMES/dataBase/lineParameter.db')
        cursor = connection.cursor()
        cursor.execute("SELECT * FROM Lines WHERE Name = ?",(name,))
        self.params =  cursor.fetchall()[0]
        self.name = self.params[1]
        self.frequency = self.params[2]
        self.wavelenght = self.params[3]

