# dao2caom2
Application to generate CAOM2 Observations from DAO FITS files.

# How To Run dao2caom2

These are Linux-centric instructions.

In an empty directory (the 'working directory'), on a machine with Docker installed:

1. Set up credentials. Due to the current state of CADC hosts, proxy certificates will not work, therefore you must use a netrc file.

    In the 'working directory', create a file named 'netrc'. 
This is the expected netrc file that will have the credentials required for the 
CADC services. These credentials allow the user to read, write, and delete 
CAOM2 observations, read file header metadata and files 
from data, and store thumbnails and previews to CADC storage. This file should have content that looks like the following:

   ```
   machine www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca login canfarusername password canfarpassword
   machine www.canfar.net login canfarusername password canfarpassword
   machine sc2.canfar.net login canfarusername password canfarpassword
   machine ws-cadc.canfar.net login canfarusername password canfarpassword
   ```
   
   1. Replace canfarusername and canfarpassword with your CADC username and 
   password values.

   1. The permissions for this file must be 600 (owner rw only):

      ```
      chmod 600 netrc
      ```
   
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

# How to Run dao2caom2 in a cron job:

In an empty directory (the 'working directory'), on a machine with Docker installed:

1. In the master branch of this repository, find the scripts directory, and copy the file dao_run_state.sh to the working directory. e.g.:

  ```
  wget https://raw.github.com/opencadc-metadata-curation/dao2caom2/master/scripts/dao_run_state.sh
  ```

2. Ensure the script is executable:

```
chmod +x dao_run_state.sh
```

3. To run the application:

```
./dao_run_state.sh
```

Note that the e-transfer script daoFileIngest knows whether or not the files should be compressed before being stored at CADC.


# Test Files

Current active instruments/detectors for "DAO" collection.  
I'm ignoring "DAOPLATES" collection which is currently static.

('c' in file name indicated a ccd is the detector; see later text...)

1. 1.2m McKellar Spectrograph

    1. dao_c122_2017_011124.fits	science observation
    1. dao_c122_2017_011118.fits	comparison arc
    1. dao_c122_2017_011111.fits	bias
    1. dao_c122_2017_011095.fits	flat
    1. dao_c122_2017_007425.fits	dark

1. 1.8m Cassegrain Spectrograph

    1. dao_c182_2017_010923.fits	science observation
    1. dao_c182_2017_010924.fits	comparison arc
    1. dao_c182_2017_010870.fits	bias
    1. dao_c182_2017_010896.fits	flat
    1. dao_c182_2017_007397.fits	dark

1. 1.8m Cassegrain Spectropolarimeter

    1. dao_c182_2017_019323.fits	science observation
    1. dao_c182_2017_019322.fits	comparison arc
    1. dao_c182_2017_019238.fits	bias
    1. dao_c182_2017_019030.fits	flat
    1. dao_c182_2017_003051.fits	dark

1. 1.8m Newtonian Imager

    1. dao_c182_2017_016292.fits	science observation
    1. dao_c182_2017_016319.fits	flat
    1. dao_c182_2017_016254.fits	bias
    1. dao_c182_2017_009814.fits	dark


-------------------------


Old data with the Reticon detector ('r' in file name).  
Limited files so far. Original headers very limited.

1. 1.2m spectra

    1. dao_r122_1989_003112.fits	science object
    1. dao_r122_1989_003115.fits	flat
    1. dao_r122_1989_003111.fits	comparison arc
    1. dao_r122_1989_003184.fits	dark

1. 1.8m spectra

    1. dao_r182_1989_000369.fits	science object
    1. dao_r182_1989_000372.fits	comparison arc
    1. dao_r182_1989_000377.fits	dark
    1. dao_r182_1989_000481.fits	flat


-------------------------


Processed data (limited 'experimental' amounts;
any suffix after odometer number indicates it's a processed file)

1. Spectra

    1. dao_c122_2007_000882_v.fits	processed science observation (suffix '_v')
        1. dao_c122_2007_000882.fits	corresponding unprocessed observation
    1. dao_c122_2007_000881_e.fits	processed arc (suffix '_e')
        1. dao_c122_2007_000881.fits	corresponding unprocessed arc
    1. dao_c122_2007_000916_F.fits	co-added, processed flat field (suffix '_F')
        1. member files are dao_c122_2007_000916.fits through dao_c122_2007_000926.fits
    1. dao_c122_2016_012652_B.fits	co-added, processed bias (suffix '_B')
        1. member files are dao_c122_2016_012652.fits to 012666 and dao_c122_2016_012726 to 012741.fits

1. Images

    1. dao_c182_2016_004034_a.fits	processed science observation (suffix '_a')
        1. dao_c182_2016_004034.fits	corresponding unprocessed observation
    1. dao_c182_2016_008332_B.fits	co-added processed bias (suffix '_B')
        1. member files dao_c182_2016_008332.fits through dao_c182_2016_008347.fits
    1. dao_c182_2016_002019_F.fits	co-added processed flat (suffix '_F')
        1. member files dao_c182_2016_002019 through dao_c182_2016_002031.fits
=======
An application to generate CAOM2 Observations from DAO FITS files.

