import gzip
# import pickle
# import numpy as np
# from sklearn.neighbors import NearestNeighbors

from ordered_set import OrderedSet

# import bz2
# import pickle 
# import _pickle as cPickle


def load_barcode_allowlist(filename):
    barcodes = OrderedSet()

    if filename is not None:
        f = gzip.open(filename, 'rb') if filename.endswith(".gz") else open(filename, 'rb')
        for l in f:
            bc = l.decode("utf-8").strip()

            barcodes.add(bc)

    return barcodes