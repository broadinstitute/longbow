import pytest

import pysam
import tempfile

from longbow.utils import bam_utils
from longbow.utils import model


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


@pytest.fixture(scope="module", params=list(model.LibraryModel.pre_configured_models.keys()))
def bam_header_with_program_group(request):
    rgid = '01234567'
    movie_name = 'm00001e_210000_000000'
    version = '0.0.0'

    model_name = request.param
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


@pytest.fixture(scope="module", params=list(filter(lambda x: x != "mas15v2", model.LibraryModel.pre_configured_models.keys())))
def bam_header_with_multiple_program_groups(request):
    rgid = '01234567'
    movie_name = 'm00001e_210000_000000'
    version = '0.0.0'

    header_dict = {
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
    }

    for model_name in ['mas15v2', request.param]:
        model_json = model.LibraryModel.build_pre_configured_model(model_name).to_json(indent=None)

        header_dict['PG'].append({
            'ID':f'longbow-annotate-{version}',
            'PN':'longbow',
            'VN':version,
            'DS':f'Annotate reads in a BAM file with segments from the model.  MODEL(s): {model_json}',
            'CL':f'longbow annotate -m {model_name} -o test.ann.bam test.bam'
        })

    header = pysam.AlignmentHeader.from_dict(header_dict)

    return header


def test_load_read_count():
    pbi_file = 'tests/test_data/mas15_test_input.bam.pbi'
    n_reads_obs = bam_utils.load_read_count(pbi_file)
    n_reads_exp = 8

    assert n_reads_obs == n_reads_exp


def test_compute_shard_offsets():
    pbi_file = 'tests/test_data/mas15_test_input.bam.pbi'

    shard_offsets_exp = {
        2: [33947648, 921501696, 921547252],
        3: [33947648,  33995961, 921533666, 921547252],
        4: [33947648,  33980810, 921501696, 921533666, 921547252],
        5: [33947648,  33980810, 921501696, 921533666, 921547252]
    }

    for num_shards in range(2, 6):
        shard_offsets_obs, _, _, _, _ = bam_utils.compute_shard_offsets(pbi_file, num_shards)

        for i in range(len(shard_offsets_exp[num_shards])):
            assert shard_offsets_obs[i] == shard_offsets_exp[num_shards][i]


def test_bam_header_has_model(bam_header_with_program_group):
    ret = bam_utils.bam_header_has_model(bam_header_with_program_group)
    assert ret == True


def test_bam_header_missing_model(bam_header_without_program_group):
    ret = bam_utils.bam_header_has_model(bam_header_without_program_group)
    assert ret == False


@pytest.mark.parametrize("bam_header_with_multiple_program_groups", ['mas15teloprimev2'], indirect=True)
def test_load_models_from_bam_header(bam_header_with_multiple_program_groups):
    with tempfile.NamedTemporaryFile(delete=True) as f:
        with pysam.AlignmentFile(f.name, "wb", header=bam_header_with_multiple_program_groups) as bf:
            bf.close()

        models = bam_utils.load_models(None, input_bam=f)

        assert len(models) == 2
        assert models[0].name == 'mas15v2'
        assert models[1].name == 'mas15teloprimev2'


def test_load_model_from_name():
    for model_name in list(model.LibraryModel.pre_configured_models.keys()):
        lb_models = bam_utils.load_models([model_name])

        assert lb_models[0].name == model_name


def test_load_model_from_json():
    with tempfile.NamedTemporaryFile(delete=False) as f:
        f.write(b'{"name":"mas3teloprimev2","description":"An example model for a 3-element array.","version":"1.0.0","array_element_structure":[["A","TPV2_adapter","cDNA","Poly_A","idx","rev_bind"],["B","TPV2_adapter","cDNA","Poly_A","idx","rev_bind"], ["C","TPV2_adapter","cDNA","Poly_A","idx","rev_bind","D"]],"adapters":{"TPV2_adapter":"CTACACGACGCTCTTCCGATCTTGGATTGATATGTAATACGACTCACTATAG","rev_bind":"CTCTGCGTTGATACCACTGCTT","A":"AGCTTACTTGTGAAGAT","B":"ACTTGTAAGCTGTCTAT","C":"ACTCTGTCAGGTCCGAT","D":"ACCTCCTCCTCCAGAAT","Poly_A":{"HomopolymerRepeat":["A",30]},"idx":{"FixedLengthRandomBases":10},"cDNA":"random"},"direct_connections":{"A":["TPV2_adapter"],"B":["TPV2_adapter"],"C":["TPV2_adapter"],"D":["TPV2_adapter"],"TPV2_adapter":["cDNA"],"cDNA":["Poly_A"],"Poly_A":["idx"],"idx":["rev_bind"],"rev_bind":["B","C","D"]},"start_element_names":["A","TPV2_adapter"],"end_element_names":["rev_bind","D"],"named_random_segments":["idx","cDNA"],"coding_region":"cDNA","annotation_segments":{"idx":[["BC","XB"]]}}\n')

    lb_models = bam_utils.load_models([f.name])

    assert lb_models[0].name == "mas3teloprimev2"


def test_reverse_complement():
    d = {
        'GTTCTAGCGCTAGTATG': 'CATACTAGCGCTAGAAC',
        'CTGCATAAAT': 'ATTTATGCAG',
        'TTTCGCGCATATAG': 'CTATATGCGCGAAA',
        'CCAA': 'TTGG',
        'TTATATTTAT': 'ATAAATATAA',
    }

    for k, v_exp in d.items():
        v_obs = bam_utils.reverse_complement(k)

        assert v_obs == v_exp


@pytest.mark.slow
def test_generate_read_name_produces_no_collisions():
    movie_name = 'm64020e_210000_000000'

    read_names = set()
    for zmw in range(8000000):
        for split_read_index in range(20):
            new_read_name = bam_utils.generate_read_name(movie_name, zmw, split_read_index)

            assert new_read_name not in read_names

            read_names.add(new_read_name)
