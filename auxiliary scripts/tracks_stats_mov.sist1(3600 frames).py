#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Count tracks & videos in excel files and supp. material
import os
from collections import defaultdict
from datetime import datetime
import time
# import re
import sys
import pandas as pd
import logging
from PIL import Image


logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(message)s",
    datefmt='%d-%b-%y %H:%M:%S',
    handlers=[
        logging.FileHandler(fr'C:\Users\pc\Desktop\Doctorado\Publicaciones\mov.sist\valid_tracks_707\Track_info_{datetime.now().strftime("%d.%m.%Y_%H.%M.%S")}.txt'),
        logging.StreamHandler(sys.stdout)
    ])


t = time.time()


walk_dir = r'C:\Users\pc\Desktop\Doctorado\Publicaciones\mov.sist\valid_tracks_707'


logging.info('\n\n\nwalk_dir = ' + walk_dir)
# If your current working directory may change during script execution, it's recommended to
# immediately convert program arguments to an absolute path. Then the variable root below will
# be an absolute path as well.
walk_dir = os.path.abspath(walk_dir)
logging.info('walk_dir (absolute) = ' + os.path.abspath(walk_dir))


excelfiles = []
tracks_quantity = defaultdict()
numreps = defaultdict()
numcells = defaultdict(defaultdict().copy)
max_min_cells = defaultdict(defaultdict().copy)
left_n = defaultdict(defaultdict().copy)
right_n = defaultdict(defaultdict().copy)
x_displacement_serie = defaultdict(defaultdict().copy)


for root, subdirs, files in os.walk(walk_dir):
    logging.info(f'\n\nNOW IN ROOT = {root}\n')
    logging.info(f'\n\n\tThere are {len(files)} files in {root}\n')
    

    for subdir in subdirs:
        logging.info('\t- subdirectory ' + subdir)


    for filename in files:
        if filename.endswith('.xlsx'):
            excelfiles.append(filename)
            scenario = root.rsplit('\\',2)[1]
            cell_type = root.rsplit('\\',1)[1]
            experiment = scenario + '_' + cell_type
            file_path = os.path.join(root, filename)
            serie = pd.read_excel(file_path, header=None)
            group_ids = (serie[7] < (serie[7].shift())).cumsum()            
            #Coordenadas X iniciales de cada célula en esa serie
            o = serie.groupby(group_ids, as_index=False).first()[4] 
            #Coordenadas X finales de cada célula en esa serie
            f = serie.groupby(group_ids, as_index=False).last()[4]
            x_displacement_serie[experiment][filename] = [f_i-o_i for o_i,f_i in zip(o.tolist(), f.tolist())]
            if any(ele == 0 for ele in x_displacement_serie[experiment][filename]): # detect cells with no displacement
                logging.critical(f'\t- cell(s) nº {[idx for idx, val in enumerate(x_displacement_serie[experiment][filename]) if val == 0]} have no displacement')
                raise SystemExit(f'\t- cell(s) nº {[idx for idx, val in enumerate(x_displacement_serie[experiment][filename]) if val == 0]} have no displacement')
            left_n[experiment][filename] = sum(1 for value in x_displacement_serie[experiment][filename] if value < 0 )
            right_n[experiment][filename] = sum(1 for value in x_displacement_serie[experiment][filename] if value > 0 )
            numcells[experiment][filename] = len(x_displacement_serie[experiment][filename])
            logging.info(f'\n\t- video {filename} contains {numcells[experiment][filename]} tracks and {len(serie)} data points,{chr(10)}{chr(10)}\t of which {left_n[experiment][filename]} end on the left side and {right_n[experiment][filename]} on the right side, with ratios {(left_n[experiment][filename]/len(x_displacement_serie[experiment][filename]))*100}% left and {(right_n[experiment][filename]/len(x_displacement_serie[experiment][filename]))*100}% right{chr(10)}{chr(10)}')


    if any(file.endswith('.xlsx') for file in files): # If the folder that was just parsed by the previous for loop contained any xlsx files do...
        max_min_cells[experiment]['max'] = max(numcells[experiment].values())
        max_min_cells[experiment]['min'] = min(numcells[experiment].values())
        logging.info(f'Max and min number of cells per rep in {experiment} are {max_min_cells[experiment]["min"]}-{max_min_cells[experiment]["max"]}')
        numreps[experiment] = len(files)
        logging.info(f'{experiment} contains {numreps[experiment]} videos')
        tracks_quantity[experiment] = sum(numcells[experiment].values())
        logging.info(f'{experiment} contains {tracks_quantity[experiment]} tracks')
        logging.info(f'Tracks in {experiment} displacement values:{chr(10)}{chr(10).join(str(v) for v in sum(x_displacement_serie[experiment].values(),[]))}')
        logging.info(f'{experiment} has {sum(right_n[experiment].values())} cells to the right')
        logging.info(f'{experiment} has {sum(left_n[experiment].values())} cells to the left')
        logging.info(f'{experiment} has {(sum(left_n[experiment].values())/tracks_quantity[experiment])*100} left side ratio')
        logging.info(f'{experiment} has {(sum(right_n[experiment].values())/tracks_quantity[experiment])*100} right side ratio')
    else:
        logging.info(f'\nThere are no xlsx files in folder: {root}')


logging.info(f'All subfolders contain {sum(numreps.values())} videos and {sum(tracks_quantity.values())} tracks')
elapsed = time.time() - t
logging.info(f'{elapsed} seconds elapsed')


# Count "leaf" folders (videos) in the supp. material directory

supp_material = r'C:\Users\pc\Desktop\Doctorado\Publicaciones\mov.sist\Supp_mat con vídeos mov sist'

logging.info(f'{chr(10)}{chr(10)}Now counting videos in Supp. Material and xlsx files in \
{supp_material}{chr(10)}{chr(10)}')
folders = []
for root, dirs, files in os.walk(supp_material):
    for file in files:
        if file.endswith(".jpg"):
            folders.append(os.path.join(root, file))
            img = Image.open(os.path.join(root, file))
            wid, hgt = img.size
            logging.info(f'{chr(10)}File {root.split(chr(92),6)[-1]} has {str(wid) + "x" + str(hgt)} resolution')
            break
logging.info(f'There are {len(folders)} videos in "{supp_material}" and {len(excelfiles)} xlsx files in "{walk_dir}"{chr(10)}')