### using script to convert svs to png or another image type

import os
import glob
from wand.image import Image
#import functools
import re
import gc

### quick implementation of paste0  from R on lists
#paste0 = functools.partial(paste, sep="")

wd = '/home/kevhu/data/TCGA-BRCA'
os.chdir(wd)



for file in glob.glob('./*.svs'):
    with Image(filename=file) as img:
        img.format = 'png'
        name = str(file).replace('svs', 'png')
        img.save(filename= name)
	gc.collect()
