# dao2caom2
An application to generate CAOM2 Observations from DAO FITS files.

# How To Run dao2caom2

These are Linux-centric instructions.

In an empty directory (the 'working directory'), on a machine with Docker installed:

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

