import pytest
from longbow.utils import bam_utils
from longbow.utils import model

import pysam


# fp = tempfile.TemporaryFile()
# with pysam.AlignmentFile(fp, mode="wb", header=header, check_header=True, check_sq=False) as bam_file:
#     read = pysam.AlignedSegment(header)
#     read.query_name = f'{movie_name}/{zmw}/ccs'
#     read.mapping_quality = 255
#     read.query_sequence  = 'AAACCCGGGTTT'
#     read.query_qualities = [10] * len(read.query_sequence)
#     read.set_tag('RG', rgid)
#     read.set_tag('np', 5)
#     read.set_tag('rq', 0.999)
#     read.set_tag('zm', zmw)
#     read.set_tag('YN', model)
#     # SG:Z:TPV2_adapter:0-58,cDNA:59-292
#     # YS:f:-462.004
#     # RC:i:1
#     # XQ:Z:20/104,0/0
#     # YQ:Z:0.1923
#     bam_file.write(read)


@pytest.fixture
def bam_header_without_program_group():
    rgid = '01234567'
    movie_name = 'm00001e_210000_000000'

    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.5', 'SO': 'unknown', 'pb': '3.0.1'},
        'RG': [{
            'ID':rgid,
            'PL':'PACBIO',
            'DS':'READTYPE=CCS;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-894-200;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000',
            'LB':'TestLib',
            'PU':movie_name,
            'SM':'TestSample',
            'PM':'SEQUELII',
            'CM':'S/P5-C2/5.0-8M'
        }],
        'PG': [{
            'ID':'ccs-6.0.0',
            'PN':'ccs',
            'VN':'6.0.0',
            'DS':'Generate circular consensus sequences (ccs) from subreads.',
            'CL':'ccs /opt/pacbio/pa-ccs/current/bin/ccs --all --streamed /data/pa/m00001e_210000_000000.consensusreadset.xml --bam /data/pa/m00001e_210000_000000.reads.bam --suppress-reports --num-threads 232 --log-level INFO --log-file /data/pa/m64020e_210000_000000.ccs.log --report-json /data/pa/m64020e_210000_000000.ccs_reports.json --report-file /data/pa/m64020e_210000_000000.ccs_reports.txt --metrics-json /data/pa/m64020e_210000_000000.zmw_metrics.json.gz --hifi-summary-json /data/pa/m64020e_210000_000000.hifi_summary.json --stderr-json-log --all-kinetics --subread-fallback'
        }],
        'SQ': []
    })

    return header


@pytest.fixture
def bam_header_with_program_group():
    rgid = '01234567'
    movie_name = 'm00001e_210000_000000'
    version = '0.0.0'

    model_name = 'mas15v2'
    model_json = model.LibraryModel.build_pre_configured_model(model_name).to_json(indent=None)

    header = pysam.AlignmentHeader.from_dict({
        'HD': {'VN': '1.5', 'SO': 'unknown', 'pb': '3.0.1'},
        'RG': [{
            'ID':rgid,
            'PL':'PACBIO',
            'DS':'READTYPE=CCS;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-894-200;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000',
            'LB':'TestLib',
            'PU':movie_name,
            'SM':'TestSample',
            'PM':'SEQUELII',
            'CM':'S/P5-C2/5.0-8M'
        }],
        'PG': [{
            'ID':'ccs-6.0.0',
            'PN':'ccs',
            'VN':'6.0.0',
            'DS':'Generate circular consensus sequences (ccs) from subreads.',
            'CL':'ccs /opt/pacbio/pa-ccs/current/bin/ccs --all --streamed /data/pa/m00001e_210000_000000.consensusreadset.xml --bam /data/pa/m00001e_210000_000000.reads.bam --suppress-reports --num-threads 232 --log-level INFO --log-file /data/pa/m64020e_210000_000000.ccs.log --report-json /data/pa/m64020e_210000_000000.ccs_reports.json --report-file /data/pa/m64020e_210000_000000.ccs_reports.txt --metrics-json /data/pa/m64020e_210000_000000.zmw_metrics.json.gz --hifi-summary-json /data/pa/m64020e_210000_000000.hifi_summary.json --stderr-json-log --all-kinetics --subread-fallback'
        },
        {
            'ID':f'longbow-annotate-{version}',
            'PN':'longbow',
            'VN':version,
            'DS':f'Annotate reads in a BAM file with segments from the model.  MODEL(s): {model_json}',
            'CL':f'longbow annotate -m {model_name} -o test.ann.bam test.bam'
        }],
        'SQ': []
    })

    return header


def test_bam_header_has_model(bam_header_with_program_group):
    ret = bam_utils.bam_header_has_model(bam_header_with_program_group)
    assert ret == True


def test_bam_header_missing_model(bam_header_without_program_group):
    ret = bam_utils.bam_header_has_model(bam_header_without_program_group)
    assert ret == False


