# dao2caom2
Application to generate CAOM2 Observations from DAO FITS files.


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
