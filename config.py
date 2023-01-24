import os

absolute = os.path.abspath(os.getcwd())
absolute = absolute + '/'

data = {'data': absolute + 'data',
        'plots': absolute + 'plots'
}
# Define Mongo configuration
db_host = '172.21.2.89'
db_port = 27017
db_name = 'biodbs'
db_coll = 'seqres'
