#!/bin/bash

# This just checks if Python3  is available, installs it if it is not,
# then calls InstallPrerequisites-Ubuntu-New.py

# Check if Python 3 is installed.
which python3 > /dev/null

# If Python 3 is not installed, install it now.
if [ $? -ne 0 ] ; then
    sudo apt install python3
fi

scriptsDirectory="$(dirname "${BASH_SOURCE}")"
${scriptsDirectory}/InstallPrerequisites-Ubuntu-New.py


