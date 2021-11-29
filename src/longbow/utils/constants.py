################################################################################
# General / high-level constants:
#######################################

DEFAULT_MODEL = "mas15v2"


################################################################################
# Constants for model construction:
#######################################

RANDOM_SEGMENT_NAME = "random"
FIXED_LENGTH_RANDOM_SEGMENT_TYPE_NAME = "FixedLengthRandomBases"
HPR_SEGMENT_TYPE_NAME = "HomopolymerRepeat"

RANDOM_SILENT_STATE_A = "RDA"
RANDOM_SILENT_STATE_B = "RDB"
RANDOM_BASE_STATE = "RI"

START_STATE_INDICATOR = "-start"
END_STATE_INDICATOR = "-end"

BAKE_MERGE_STRATEGY = "None"
MAS_SCAFFOLD_NAMES = {"MARS", "VENUS", "ZEPHYR", "AUSTER", "BOREAS"}


################################################################################
# Constants for bam file reading / writing:
#######################################

SEGMENTS_TAG = "SG"
SEGMENTS_QUAL_TAG = "XQ"
SEGMENTS_RC_TAG = "RC"
SEGMENT_TAG_DELIMITER = ","

READ_IS_SEGMENTED_TAG = "ZS"
READ_MODEL_NAME_TAG = "YN"
READ_MODEL_SCORE_TAG = "YS"
READ_IS_VALID_FOR_MODEL_TAG = "YV"
READ_FIRST_KEY_SEG_TAG = "YK"
READ_NUM_KEY_SEGMENTS_TAG = "YG"
READ_APPROX_QUAL_TAG = "YQ"
READ_ADAPTER_TAG = 'ZA'
READ_ADAPTER_POS_TAG = "XA"
READ_UMI_TAG = 'ZU'
READ_UMI_POS_TAG = "XU"
READ_RAW_UMI_TAG = "XM"  # UMI sequence (for IsoSeq3 compatibility - https://isoseq.how/general-faq.html)
READ_BARCODE_TAG = 'CR'  # Cell barcode
READ_RAW_BARCODE_TAG = "XC"  # barcode sequence (for IsoSeq3 compatibility - https://isoseq.how/general-faq.html)
READ_BARCODE_POS_TAG = "XB"
READ_BARCODE_QUAL_TAG = "CY"  # Cell barcode read quality
READ_BARCODE_CORRECTED_TAG = 'CB'  # Cell barcode that is error-corrected and confirmed against a list of known-good barcode sequences
READ_BARCODE_CONF_FACTOR_TAG = "XF"
READ_TAGS_ORDER_TAG = "XA"  # Order of tag names
READ_NUM_CONSENSUS_PASSES_TAG = "ic"  # Sum of number of passes from all ZMWs used to create consensus (e.g. 1)
READ_ZMW_NAMES_TAG = "im"  # ZMW names associated with this isoform (e.g. m64013e_211031_055434/1/ccs)
READ_NUM_ZMWS_TAG = "is"  # Number of ZMWs associated with this isoform (e.g. 1)
READ_CLIPPED_SEQS_LIST_TAG = "it"  # List of barcodes/UMIs clipped during tag (e.g. TCAGGTGCAGGTCGGATCCTGCGCAT)
READ_ZMW_TAG = "zm"
READ_ALTERED_NAME_TAG = "XN"  # Altered read name given by Longbow to a segmented read (used for debugging)
CONF_FACTOR_SCALE = 100
READ_SPATIAL_BARCODE1_TAG = "X1"
READ_SPATIAL_BARCODE1_POS_TAG = "XP"
READ_SPATIAL_BARCODE2_TAG = "X2"
READ_SPATIAL_BARCODE2_POS_TAG = "XQ"