#! /bin/bash
/bin/rm -rf /www/GeoCLUSTER.orig
/bin/mv /www/GeoCLUSTER /www/GeoCLUSTER.orig
/bin/mv /www/GeoCLUSTER.new /www/GeoCLUSTER
systemctl restart httpd