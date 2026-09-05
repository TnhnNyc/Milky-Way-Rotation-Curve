#!/usr/bin/env python
# coding: utf-8

# # Query
# ---
# Functions to fetch **Gaia DR3** astrometry and **APOGEE DR17/SDSS** spectroscopic
# parameters.
# 
# NOTE: These functions need internet access to the Gaia archive (gea.esac.esa.int) 
# and the SDSS Science Archive Server (data.sdss.org).
# 

from __future__ import annotations
import os

def fetch_apogee_dataset(dr=17, local_dir="data"):
    
    # ------------------------------------------------------------
    # 1. Prepare the local cache path
    # ------------------------------------------------------------
    filename = f"apogee_astroNN-DR{dr}.fits"
    os.makedirs(local_dir, exist_ok=True)
    destination = os.path.join(local_dir, filename)

    # ------------------------------------------------------------
    # 2. Data file control
    # ------------------------------------------------------------
    if os.path.exists(destination):
        print(f'APOGEE DR{dr} data is already found: {destination}')
        return load_apogee_table(destination)
    else:
        print('Downloading APOGEE data will start...')

    # ------------------------------------------------------------
    # 3. Querying APOGEE data
    # ------------------------------------------------------------  
    
    # downloading the data via astroNN library (recommended)
    try:
        from astroNN.apogee import apogee_astronn

        print(f"Downloading APOGEE DR{dr} data via astroNN...")
        data_table = apogee_astronn(dr=dr) 

        data_table.write(destination, format='fits', overwrite=True)
        print(f"Data successfully cached to {destination}")
        return data_table

    # a direct download from the SDSS Science Archive Server (it runs if astroNN library is not found)
    except ImportError:
        import urllib.request
        
        url = f"https://data.sdss.org/sas/dr{dr}/apogee/vac/apogee-astronn/{filename}"
        print(f"astroNN package not found; downloading directly from:\n  {url}")
        urllib.request.urlretrieve(url, destination)
        print(f"Data successfully downloaded and saved to {destination}")
        return load_apogee_table(destination)

def load_apogee_table(fits_path):
    """Load the astroNN VAC FITS table into an astropy Table."""
    from astropy.table import Table
    return Table.read(fits_path)

def fetch_gaia_dr3_dataset(source_ids=None, batch_size=3000):
    from astroquery.gaia import Gaia
    from astropy.table import vstack, Table
    import math

    # ------------------------------------------------------------
    # 1. Prepare the local cache path
    # ------------------------------------------------------------
    os.makedirs('data', exist_ok=True)
    output_path = os.path.join('data', 'gaia_dr3.fits')

    temp_dir = os.path.join('data', 'temp_batches')
    os.makedirs(temp_dir, exist_ok=True)

    # ------------------------------------------------------------
    # 2. Data file control
    # ------------------------------------------------------------
    if os.path.exists(output_path):
        print(f'Gaia DR3 data is already found: {output_path}')
        return Table.read(output_path, format='fits')
    else:
        print('Downloading Gaia DR3...')
        
    # ------------------------------------------------------------
    # 3. Querying Gaia DR3 data
    # ------------------------------------------------------------    
    if source_ids is not None: 
        
        # Counter of batch
        total_batches = math.ceil(len(source_ids)/batch_size)
        current_batch = 0
            
        for start in range(0, len(source_ids), batch_size):
            # Counter:
            current_batch += 1
                        
            # Preparing a file path to save current batch and control it whether exist or not before 
            temp_batch_path = os.path.join(temp_dir, f"batch_{current_batch:04d}.fits")

            if os.path.exists(temp_batch_path):
                print(f"Batch {current_batch}/{total_batches} already exists locally. Skipping...", end="\r")
                continue
            
            # Processing: 
            print(f'Processing batch: {current_batch}/{total_batches}...', end='\r')
            batch = source_ids[start:start + batch_size]

            # the source_ids elements have an int structure, so they should be converted string structure
            # this process is neccessary to create an ADQL query to fetch the Gaia DR3 dataset:
            ids_str = ",".join(str(int(s)) for s in batch)

            # creating an ADQL query to fetch the Gaia DR3 dataset:
            query = f"""
            SELECT
                source_id,
                ra, dec,
                parallax, parallax_error,
                pmra, pmra_error,
                pmdec, pmdec_error,
                radial_velocity, radial_velocity_error,
                phot_g_mean_mag, ruwe
            FROM gaiadr3.gaia_source
            WHERE source_id IN ({ids_str})
            """
            # NOTE: We only get some properties which are written in the SELECT part of ADQL query

            # Submitting the ADQL guert to the Gaia Archive Server:
            try:
                job = Gaia.launch_job_async(query)
                result = job.get_results()
                # Saving the current batch to the temp_batches file
                result.write(temp_batch_path, format='fits', overwrite=True)
            except Exception as e:
                print(f"\n[ERROR] Connection lost or server issue at batch {current_batch}!")
                print("Your progress has been saved. Please run the code again later to resume.")
                raise e 

        # Merging all the batch files
        print("\nAll batches downloaded successfully! Merging temporary files...")

        temp_files = sorted([os.path.join(temp_dir, i) for i in os.listdir(temp_dir) if i.endswith('.fits')])
        tables = []
        for file_path in temp_files:
            tables.append(Table.read(file_path, format='fits'))

        final_table = vstack(tables)

        # Saving Gaia DR3 dataset:
        final_table.write(output_path, format='fits', overwrite=True)
        print(f'Data successfully saved to {output_path}')

        """
        # Cleaning temporary .fits files in the temp_batches folder
        for file_path in temp_files:
            os.remove(file_path)
        os.rmdir(temp_dir)
        print("Temporary batch files cleaned up.")
        """
        return final_table
    else:
        raise ValueError('Provide either surce_ids or ra_deg/dec_deg/radius_deg')

    
