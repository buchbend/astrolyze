import socket
computerName = socket.gethostname()

if computerName == 'lt-cb':
    gralx1 = '/gra-lx1/users/buchbend/'
    sys = '/home/buchbend/'
    prog = '/home/buchbend/Programs/'
    wcs = '/home/buchbend/Programs/wcstools-3.8.4/bin/'
    dataBase = '/home/buchbend/HERMES/dataBase/'

if computerName == 'gra-lx18':
    gralx1 = '/gra-lx1/users/buchbend/'
    sys = '/home/buchbend/'
    prog = '/home/buchbend/Programs/'
    wcs = '/home/buchbend/Programs/wcstools-3.8.4/bin/' 
    dataBase = '/home/buchbend/HERMES/dataBase/'
