#! /bin/bash
/bin/rm -rf /www/GeoCLUSTER.orig

# activating the virtual environment
/www/geoVenv/bin/pip install --upgrade pip
/www/geoVenv/bin/pip install -r /www/GeoCLUSTER.new/dash_app/requirements.txt

# downloading hdf5 if it is not already downloaded
# cd /www/GeoCLUSTER.new/dash_app/data; /www/geoVenv/bin/python download_hdf5.py

#giving apache permission to the venv
chown -R apache:apache /www/geoVenv


/bin/mv /www/GeoCLUSTER /www/GeoCLUSTER.orig
/bin/mv /www/GeoCLUSTER.new /www/GeoCLUSTER
echo "deployment_type=aws" > /www/GeoCLUSTER/dash_app/.env
chown -Rh apache /www/GeoCLUSTER/dash_app
chmod -R o+rx /www/GeoCLUSTER
systemctl restart httpd