import os
USER = os.getenv("USER")
print USER

prefifx_database = "/home/{}/.astrolyze/database/".format(USER)
