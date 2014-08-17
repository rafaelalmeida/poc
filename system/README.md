SETUP INSTRUCTIONS

Put the 2014 full set contest data inside the /system/data folder.

DEPENDENCIES

- GDAL (with Python bindings)
- OpenCV

INSTALLING GDAL

Cygwin

Build from source. Download source, then

`./configure --with-libtool=no --with-python`

The --with-libtool=no option makes compilation a little faster on 
Cygwin.

Linux (tested on Ubuntu)

`sudo apt-get install libgdal`
`sudo apt-get install python-gdal`
