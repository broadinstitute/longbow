import gzip
from ordered_set import OrderedSet


def load_barcode_allowlist(filename):
    barcodes = OrderedSet()

    if filename is not None:
        f = gzip.open(filename, 'rb') if filename.endswith(".gz") else open(filename, 'rb')
        for l in f:
            bc = l.decode("utf-8").strip()

            barcodes.add(bc)

    return barcodes