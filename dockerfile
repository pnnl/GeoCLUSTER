# Dockerfile

# Choose a Python base image. 'slim' versions are smaller.
FROM python:3.8.10-slim

# Set environment variables to prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install Apache and the WSGI module FOR PYTHON 3.8, plus cleanup
# --- MODIFIED LINE BELOW ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    apache2 \
    apache2-dev \
    python3-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set the working directory inside the container
WORKDIR /app

#building mod wsgi
#copy
COPY mod_wsgi-5.0.2.tar.gz ./
RUN tar xvfz mod_wsgi-5.0.2.tar.gz && \
    cd mod_wsgi-5.0.2 && \ 
    ./configure --with-apxs=/usr/bin/apxs --with-python=/usr/local/bin/python && \
    make && make install

# Copy the requirements file first to leverage Docker layer caching
COPY requirements.txt ./

# Install Python dependencies
# --no-cache-dir keeps the image size down
RUN pip install --upgrade pip
RUN pip install setuptools
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of your application code into the container
COPY . .

# Copy the custom Apache configuration to replace the default site
# This assumes your config file is named 'apache-config.conf' in the build context
COPY apache-config.conf /etc/apache2/sites-available/000-default.conf

# Add this RUN command to fix permissions
RUN chown -R www-data:www-data /app && \
    find /app -type d -exec chmod 755 {} + && \
    find /app -type f -exec chmod 644 {} + && \
    chmod +x /app/dash_app/wsgi.py # Ensure wsgi script is executable if needed (though read usually suffices)
    

#downloading python file
RUN python /app/dash_app/data/download_hdf5.py


# Expose port 80 (Apache's default HTTP port)
EXPOSE 80

# Define the command to run Apache in the foreground when the container starts
# This keeps the container running
CMD ["apache2ctl", "-D", "FOREGROUND"]