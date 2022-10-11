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
class Dataset(object):
    #def __init__(self, fname, img_transform=None, mask_transform=None, edge_weight=False):
    ### no edges or masks so getttinf rid of a few things until it reflects what we have
    def __init__(self, fname):
        # nothing special here, just internalizing the constructor parameters
        self.fname = fname
        #self.edge_weight = edge_weight

        #self.img_transform = img_transform
        #self.mask_transform = mask_transform

        self.tables = tables.open_file(self.fname)
        self.numpixels = self.tables.root.numpixels[:]
        self.nitems = self.tables.root.img.shape[0]
        self.tables.close()

        self.img = None
        #self.mask = None

    def __getitem__(self, index):
        # opening should be done in __init__ but seems to be
        # an issue with multithreading so doing here
        with tables.open_file(self.fname, 'r') as db:
            self.img = db.root.img
            ### no masks so getting rid of it first
            #self.mask = db.root.mask

            # get the requested image and mask from the pytable
            img = self.img[index, :, :, :]
            #mask = self.mask[index, :, :]

        # the original Unet paper assignes increased weights to the edges of the annotated objects
        # their method is more sophistocated, but this one is faster, we simply dilate the mask and
        # highlight all the pixels which were "added"
        #if (self.edge_weight):
        #    weight = scipy.ndimage.morphology.binary_dilation(mask == 1, iterations=2) & ~mask
        #else:  # otherwise the edge weight is all ones and thus has no affect
        #    weight = np.ones(mask.shape, dtype=mask.dtype)

        #mask = mask[:, :, None].repeat(3, axis=2)  # in order to use the transformations given by torchvision
        #weight = weight[:, :, None].repeat(3,
        #                                   axis=2)  # inputs need to be 3D, so here we convert from 1d to 3d by repetition

        img_new = img
        #mask_new = mask
        #weight_new = weight

        #seed = random.randrange(sys.maxsize)  # get a random seed so that we can reproducibly do the transofrmations
        #if self.img_transform is not None:
        #    random.seed(seed)  # apply this seed to img transforms
        #    img_new = self.img_transform(img)

        #if self.mask_transform is not None:
        #    random.seed(seed)
        #    mask_new = self.mask_transform(mask)
        #    mask_new = np.asarray(mask_new)[:, :, 0].squeeze()

        #    random.seed(seed)
        #    weight_new = self.mask_transform(weight)
        #    weight_new = np.asarray(weight_new)[:, :, 0].squeeze()

        #return img_new, mask_new, weight_new
        return img_new

    def __len__(self):
        return self.nitems

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







"""
for base_fname in files:
    base_fname = base_fname.strip()
    for fname in sorted(glob.glob("%s_*.tif" % (base_fname))):  # get all of the files which start with this patient ID
        print
        fname

        newfname_class = "%s/%s_class.png" % (OUTPUT_DIR, fname[0:-4])  # create the new files
        newfname_prob = "%s/%s_prob.png" % (OUTPUT_DIR, fname[0:-4])

        if (os.path.exists(
                newfname_class)):  # if this file exists, skip it...this allows us to do work in parallel across many machines
            continue
        print
        "working on file: \t %s" % fname

        outputimage = np.zeros(shape=(10, 10))
        scipy.misc.imsave(newfname_class,
                          outputimage)  # first thing we do is save a file to let potential other workers know that this file is being worked on and it should be skipped

        image = caffe.io.load_image(fname)  # load our image
        #		image = caffe.io.resize_image(image, [image.shape[0]/2,image.shape[1]/2]) #if you need to resize or crop, do it here
        image = np.lib.pad(image, ((hwsize, hwsize), (hwsize, hwsize), (0, 0)),
                           'symmetric')  # mirror the edges so that we can compute the full image

        outputimage_probs = np.zeros(
            shape=(image.shape[0], image.shape[1], 3))  # make the output files where we'll store the data
        outputimage_class = np.zeros(shape=(image.shape[0], image.shape[1]))
        for rowi in xrange(hwsize + 1, image.shape[0] - hwsize):
            print
            "%s\t (%.3f,%.3f)\t %d of %d" % (
            fname, time.time() - start_time, time.time() - start_time_iter, rowi, image.shape[0] - hwsize)
            start_time_iter = time.time()
            patches = []  # create a set of patches, oeprate on a per column basis
            for coli in xrange(hwsize + 1, image.shape[1] - hwsize):
                patches.append(image[rowi - hwsize:rowi + hwsize, coli - hwsize:coli + hwsize, :])

            prediction = net.predict(patches)  # predict the output
            pclass = prediction.argmax(axis=1)  # get the argmax
            outputimage_probs[rowi, hwsize + 1:image.shape[1] - hwsize,
            0:2] = prediction  # save the results to our output images
            outputimage_class[rowi, hwsize + 1:image.shape[1] - hwsize] = pclass

        outputimage_probs = outputimage_probs[hwsize:-hwsize, hwsize:-hwsize, :]  # remove the edge padding
        outputimage_class = outputimage_class[hwsize:-hwsize, hwsize:-hwsize]

        scipy.misc.imsave(newfname_prob, outputimage_probs)  # save the files
        scipy.misc.imsave(newfname_class, outputimage_class)
"""
