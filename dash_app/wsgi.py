import sys
import site

# sourced scripts
# from paths import absolute_path

# path to virtual environment
site.addsitedir('/opt/www/GeoCLUSTER/lib/python3.8/site-packages') # updated for openei
# path to app directory
sys.path.insert(0, "/var/www/GeoCLUSTER") # updated for openei; same as absolute_path

# from app import app as applicaton
from app import server as application # this targets the Flask server of the Dash app; aws