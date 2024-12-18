#! /bin/bash
/bin/rm -rf /www/GeoCLUSTER.orig
/bin/mv /www/GeoCLUSTER /www/GeoCLUSTER.orig
/bin/mv /www/GeoCLUSTER.new /www/GeoCLUSTER
chown -Rh apache /www/GeoCLUSTER/dash_app
chmod o+rx /www/GeoCLUSTER
systemctl restart httpd