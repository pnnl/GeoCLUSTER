#! /bin/bash
/bin/rm -rf /www/GeoCLUSTER.orig
cd /www/GeoCLUSTER.new/dash_app/data; python3.11 download_hdf5.py
/bin/mv /www/GeoCLUSTER /www/GeoCLUSTER.orig
/bin/mv /www/GeoCLUSTER.new /www/GeoCLUSTER
chown -Rh apache /www/GeoCLUSTER/dash_app
chmod -R o+rx /www/GeoCLUSTER
systemctl restart httpd