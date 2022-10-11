#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import scipy.misc
import time
import argparse
import sys

caffe_root = 'CAFFE ROOT'

sys.path.insert(0, caffe_root + 'python')
import caffe

# parse the command line arguments
parser = argparse.ArgumentParser(description='Generate full image outputs from Caffe DL models.')
parser.add_argument('base', help="Base Directory type {nuclei, epi, mitosis, etc}")

### KH: code expects a specific fold to run model on, but need to alter it to just work on the hdf5 table
#parser.add_argument('fold', type=int, help="Which fold to generate outputfor")
args = parser.parse_args()

# this window size needs to be exactly the same size as that used to extract the patches from the matlab version
### KH: janowczyk used mathlab script to create pathces - takes too many inodes and hard to use downstream ...
### so instead I used a pytable hdf5 format which does it all in memory

wsize = 32
hwsize = int(wsize / 2)

BASE = args.base
#FOLD = args.fold

# Locations of the necessary files are all assumed to be in subdirectoires of the base file
MODEL_FILE = '%s/models/deploy_train32.prototxt' % BASE
PRETRAINED = '%s/models/caffenet_train_w32_iter_600000.caffemodel' % (BASE)
IMAGE_DIR = '%s/images/' % BASE
OUTPUT_DIR = '%s/images/' % (BASE)

# if our output directory doesn't exist, lets create it. each fold gets a numbered directory inside of the image directory
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# load our mean file and reshape it accordingly
a = caffe.io.caffe_pb2.BlobProto()
file = open('%s/DB_train_w32.binaryproto' % (BASE), 'rb')
data = file.read()
a.ParseFromString(data)
means = a.data
means = np.asarray(means)
means = means.reshape(3, 32, 32)

# make sure we use teh GPU otherwise things will take a very long time
caffe.set_mode_cpu()
#caffe.set_mode_gpu()

# load the model
net = caffe.Classifier(MODEL_FILE, PRETRAINED,
                       mean=means,
                       channel_swap=(2, 1, 0),
                       raw_scale=255,
                       image_dims=(32, 32))

# see which files we need to produce output for in this fold
# we look at the parent IDs in the test file and only compute those images
# as they've been "held out" from the training set
#files = open('%s/test_w32_parent_%d.txt' % (BASE, FOLD), 'rb')

### KH: this is where I need to open the hdf5 table and iterate through it. below will include a lot of changes
### KH: one of the bigger parts is creating a class to manage the loading of the hdf5 table
### KH: this defines our dataset class which will be used by the dataloader. should test dataloader later before I run the script.



"""

start_time = time.time()
start_time_iter = 0

# go into the image directory so we can use glob a bit easier
os.chdir(IMAGE_DIR)

#dataset={}
###  pretty cool thing to learn about python literal "f-string", better than using the normal %
#dataset[phase]=Dataset(f"./{dataname}_{phase}.pytable")
#dataset=Dataset("./cropped_memmapSmallest_.pytable")


#output = np.zeros((0,checkpoint["n_classes"],patch_size,patch_size))
outputimage_probs = np.zeros(shape=(image.shape[0],image.shape[1],3))
outputimage_class = np.zeros(shape=(image.shape[0],image.shape[1]))


###grab the shape of the pytable
io_shape_orig = np.array(io.shape)


