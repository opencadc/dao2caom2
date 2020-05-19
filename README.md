# dao2caom2
An application to generate CAOM2 Observations from DAO FITS files.

# How To Run dao2caom2

These are Linux-centric instructions.

In an empty directory (the 'working directory'), on a machine with Docker installed:

1. Set up credentials. Due to the current state state of CADC hosts, proxy certificates will not work, therefore you must use a netrc file.

    In the 'working directory', create a file named 'netrc'. 
This is the expected netrc file that will have the credentials required for the 
CADC services. These credentials allow the user to read, write, and delete 
CAOM2 observations, and read file header metadata and files 
from data. This file should have content that looks like the following:

   ```
   machine www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca login canfarusername password canfarpassword
   machine www.canfar.net login canfarusername password canfarpassword
   machine sc2.canfar.net login canfarusername password canfarpassword
   machine ws-cadc.canfar.net login canfarusername password canfarpassword
   ```
   
   1. Replace canfarusername and canfarpassword with your CADC username and 
   password values.

   1. The permissions for this file must be 600 (owner rw only).
   
   1. The man page for netrc:
   https://www.systutorials.com/docs/linux/man/5-netrc/
   
   1. The name and location of this file may be changed by modifying the 
   netrc_filename entry in the config.yml file. This entry requires a 
   fully-qualified pathname.
1. In the master branch of this repository, one time only, find the scripts directory, and copy the file dao_run.sh to the working directory. e.g.:

   ```
   wget https://raw.github.com/opencadc-metadata-curation/dao2caom2/master/scripts/dao_run.sh
   ```

1. Ensure the script is executable, one time only:

   ```
   chmod +x dao_run.sh
   ```

1. To run the application:

    ```
    ./dao_run.sh
    ```
    
1. The config.yml file will be created in the 'working directory'. This file controls the execution of the application. See [here](https://github.com/opencadc-metadata-curation/collection2caom2/wiki/config.yml) for a description of the entries in this file.

1. To debug the application from inside the container:

   ```
   user@dockerhost:<cwd># docker run --rm -ti -v ${PWD}:/usr/src/app --name dao_run opencadc/dao2caom2 /bin/bash
   root@53bef30d8af3:/usr/src/app# dao_run
   ```

1. For some instructions that might be helpful on using containers, see:
https://github.com/opencadc-metadata-curation/collection2caom2/wiki/Docker-and-Collections

1. For some insight into what's happening, see: https://github.com/opencadc-metadata-curation/collection2caom2

