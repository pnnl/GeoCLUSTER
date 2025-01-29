import sys
import site

# sourced scripts
# from paths import absolute_path

# path to virtual environment
site.addsitedir('/var/www/geoenv-latest/lib/python3.8/site-packages') # CHANGE HERE

# path to app directory
sys.path.insert(0, "/www/GeoCLUSTER/dash_app") # CHANGE HERE; same as absolute_path

# from app import app as applicaton
from app import server as application # this targets the Flask server of the Dash app; aws