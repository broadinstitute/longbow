import itertools
import json
import pathlib
import tempfile

import pysam
import pytest

from longbow.utils import bam_utils, model, model_utils

TEST_DATA_FOLDER = pathlib.Path(__file__).parent.parent / "test_data" / "models"


@pytest.fixture
def bam_header_without_program_group():
    rgid = "01234567"
    movie_name = "m00001e_210000_000000"

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.5", "SO": "unknown", "pb": "3.0.1"},
            "RG": [
                {
                    "ID": rgid,
                    "PL": "PACBIO",
                    "DS": "READTYPE=CCS;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-894-200;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000",
                    "LB": "TestLib",
                    "PU": movie_name,
                    "SM": "TestSample",
                    "PM": "SEQUELII",
                    "CM": "S/P5-C2/5.0-8M",
                }
            ],
            "PG": [
                {
                    "ID": "ccs-6.0.0",
                    "PN": "ccs",
                    "VN": "6.0.0",
                    "DS": "Generate circular consensus sequences (ccs) from subreads.",
                    "CL": "ccs /opt/pacbio/pa-ccs/current/bin/ccs --all --streamed /data/pa/m00001e_210000_000000.consensusreadset.xml --bam /data/pa/m00001e_210000_000000.reads.bam --suppress-reports --num-threads 232 --log-level INFO --log-file /data/pa/m64020e_210000_000000.ccs.log --report-json /data/pa/m64020e_210000_000000.ccs_reports.json --report-file /data/pa/m64020e_210000_000000.ccs_reports.txt --metrics-json /data/pa/m64020e_210000_000000.zmw_metrics.json.gz --hifi-summary-json /data/pa/m64020e_210000_000000.hifi_summary.json --stderr-json-log --all-kinetics --subread-fallback",
                }
            ],
            "SQ": [],
        }
    )

    return header


@pytest.fixture(
    scope="module",
    params=list(
        map(
            lambda x: f"{x[0]}+{x[1]}",
            list(
                itertools.product(
                    list(
                        model_utils.ModelBuilder.pre_configured_models["array"].keys()
                    ),
                    list(model_utils.ModelBuilder.pre_configured_models["cdna"].keys()),
                )
            ),
        )
    ),
)
def bam_header_with_program_group(request):
    rgid = "01234567"
    movie_name = "m00001e_210000_000000"
    version = "0.0.0"

    model_name = request.param
    model_json = model.LibraryModel.build_pre_configured_model(model_name).to_json(
        indent=None
    )

    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.5", "SO": "unknown", "pb": "3.0.1"},
            "RG": [
                {
                    "ID": rgid,
                    "PL": "PACBIO",
                    "DS": "READTYPE=CCS;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-894-200;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000",
                    "LB": "TestLib",
                    "PU": movie_name,
                    "SM": "TestSample",
                    "PM": "SEQUELII",
                    "CM": "S/P5-C2/5.0-8M",
                }
            ],
            "PG": [
                {
                    "ID": "ccs-6.0.0",
                    "PN": "ccs",
                    "VN": "6.0.0",
                    "DS": "Generate circular consensus sequences (ccs) from subreads.",
                    "CL": "ccs /opt/pacbio/pa-ccs/current/bin/ccs --all --streamed /data/pa/m00001e_210000_000000.consensusreadset.xml --bam /data/pa/m00001e_210000_000000.reads.bam --suppress-reports --num-threads 232 --log-level INFO --log-file /data/pa/m64020e_210000_000000.ccs.log --report-json /data/pa/m64020e_210000_000000.ccs_reports.json --report-file /data/pa/m64020e_210000_000000.ccs_reports.txt --metrics-json /data/pa/m64020e_210000_000000.zmw_metrics.json.gz --hifi-summary-json /data/pa/m64020e_210000_000000.hifi_summary.json --stderr-json-log --all-kinetics --subread-fallback",
                },
                {
                    "ID": f"longbow-annotate-{version}",
                    "PN": "longbow",
                    "VN": version,
                    "DS": f"Annotate reads in a BAM file with segments from the model.  MODEL(s): {model_json}",
                    "CL": f"longbow annotate -m {model_name} -o test.ann.bam test.bam",
                },
            ],
            "SQ": [],
        }
    )

    return header


@pytest.fixture(
    scope="module",
    params=list(
        filter(
            lambda x: x != "mas_15+sc_10x5p",
            map(
                lambda x: f"{x[0]}+{x[1]}",
                list(
                    itertools.product(
                        list(
                            model_utils.ModelBuilder.pre_configured_models[
                                "array"
                            ].keys()
                        ),
                        list(
                            model_utils.ModelBuilder.pre_configured_models[
                                "cdna"
                            ].keys()
                        ),
                    )
                ),
            ),
        )
    ),
)
def bam_header_with_multiple_program_groups(request):
    rgid = "01234567"
    movie_name = "m00001e_210000_000000"
    version = "0.0.0"

    header_dict = {
        "HD": {"VN": "1.5", "SO": "unknown", "pb": "3.0.1"},
        "RG": [
            {
                "ID": rgid,
                "PL": "PACBIO",
                "DS": "READTYPE=CCS;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=101-894-200;SEQUENCINGKIT=101-826-100;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000",
                "LB": "TestLib",
                "PU": movie_name,
                "SM": "TestSample",
                "PM": "SEQUELII",
                "CM": "S/P5-C2/5.0-8M",
            }
        ],
        "PG": [
            {
                "ID": "ccs-6.0.0",
                "PN": "ccs",
                "VN": "6.0.0",
                "DS": "Generate circular consensus sequences (ccs) from subreads.",
                "CL": "ccs /opt/pacbio/pa-ccs/current/bin/ccs --all --streamed /data/pa/m00001e_210000_000000.consensusreadset.xml --bam /data/pa/m00001e_210000_000000.reads.bam --suppress-reports --num-threads 232 --log-level INFO --log-file /data/pa/m64020e_210000_000000.ccs.log --report-json /data/pa/m64020e_210000_000000.ccs_reports.json --report-file /data/pa/m64020e_210000_000000.ccs_reports.txt --metrics-json /data/pa/m64020e_210000_000000.zmw_metrics.json.gz --hifi-summary-json /data/pa/m64020e_210000_000000.hifi_summary.json --stderr-json-log --all-kinetics --subread-fallback",
            }
        ],
        "SQ": [],
    }

    for model_name in ["mas_15+sc_10x5p", request.param]:
        model_json = model.LibraryModel.build_pre_configured_model(model_name).to_json(
            indent=None
        )

        header_dict["PG"].append(
            {
                "ID": f"longbow-annotate-{version}",
                "PN": "longbow",
                "VN": version,
                "DS": f"Annotate reads in a BAM file with segments from the model.  MODEL(s): {model_json}",
                "CL": f"longbow annotate -m {model_name} -o test.ann.bam test.bam",
            }
        )

    header = pysam.AlignmentHeader.from_dict(header_dict)

    return header


def test_load_read_count():
    pbi_file = "tests/test_data/mas15_test_input.bam.pbi"
    n_reads_obs = bam_utils.load_read_count(pbi_file)
    n_reads_exp = 8

    assert n_reads_obs == n_reads_exp


def test_compute_shard_offsets():
    pbi_file = "tests/test_data/mas15_test_input.bam.pbi"

    shard_offsets_exp = {
        2: [33947648, 921501696, 921547252],
        3: [33947648, 33995961, 921533666, 921547252],
        4: [33947648, 33980810, 921501696, 921533666, 921547252],
        5: [33947648, 33980810, 921501696, 921533666, 921547252],
    }

    for num_shards in range(2, 6):
        shard_offsets_obs, _, _, _, _ = bam_utils.compute_shard_offsets(
            pbi_file, num_shards
        )

        for i in range(len(shard_offsets_exp[num_shards])):
            assert shard_offsets_obs[i] == shard_offsets_exp[num_shards][i]


def test_bam_header_has_model(bam_header_with_program_group):
    ret = bam_utils.bam_header_has_model(bam_header_with_program_group)
    assert ret is True


def test_bam_header_is_missing_model(bam_header_without_program_group):
    ret = bam_utils.bam_header_has_model(bam_header_without_program_group)
    assert ret is False


def _compare_models(prebuilt_model, stored_model):
    assert prebuilt_model.name == stored_model.name
    assert prebuilt_model.description == stored_model.description
    assert prebuilt_model.array_version == stored_model.array_version
    assert prebuilt_model.cdna_version == stored_model.cdna_version
    assert (
        prebuilt_model.array_element_structure == stored_model.array_element_structure
    )
    assert prebuilt_model.key_adapter_set == stored_model.key_adapter_set
    assert prebuilt_model.key_adapters == stored_model.key_adapters

    assert prebuilt_model.adapter_dict.keys() == stored_model.adapter_dict.keys()
    for p, s in zip(
        prebuilt_model.adapter_dict.values(), stored_model.adapter_dict.values()
    ):
        if isinstance(p, dict):
            assert p.keys() == s.keys()

            for k in p.keys():
                if isinstance(p[k], tuple):
                    # TODO: This is specifically to fix an issue comparing HomopolymerRepeat objects.
                    # Apparently, retrieving this from the pre-built model yields ('A', 30), but
                    # loading from json yields ['A', 30].  This should be fixed.
                    assert tuple(p[k]) == tuple(s[k])
                else:
                    assert p[k] == s[k]


def test_load_model_from_name():
    for array_model_name in list(
        model_utils.ModelBuilder.pre_configured_models["array"].keys()
    ):
        for cdna_model_name in list(
            model_utils.ModelBuilder.pre_configured_models["cdna"].keys()
        ):
            model_name = f"{array_model_name}+{cdna_model_name}"

            lb_model = bam_utils.load_model(model_name)

            stored_model = model.LibraryModel.from_json_file(
                TEST_DATA_FOLDER / f"{model_name}.json"
            )

            _compare_models(lb_model, stored_model)


def test_load_model_from_json():
    json_string = (
        b'{"name":"mas_3+bulk_teloprimeV2","description":"3-element MAS-ISO-seq array, Lexogen TeloPrime V2'
        b' kit","array":{"description":"3-element MAS-ISO-seq array","version":"3.0.0","structure":["A","B"'
        b',"C","D"],"adapters":{"A":"AGCTTACTTGTGAAGA","B":"ACTTGTAAGCTGTCTA","C":"ACTCTGTCAGGTCCGA","D":'
        b'"ACCTCCTCCTCCAGAA"},"deprecated":false,"name":"mas_3"},"cdna":{"description":"Lexogen TeloPrime V2 kit",'
        b'"version":"3.0.0","structure":["TPV2_adapter","cDNA","Poly_A","idx","rev_bind"],"adapters":'
        b'{"TPV2_adapter":"CTACACGACGCTCTTCCGATCTTGGATTGATATGTAATACGACTCACTATAG","cDNA":"random","Poly_A"'
        b':{"HomopolymerRepeat":["A",30]},"idx":{"FixedLengthRandomBases":10},"rev_bind":'
        b'"CTCTGCGTTGATACCACTGCTT"},"named_random_segments":["idx","cDNA"],"coding_region":"cDNA",'
        b'"annotation_segments":{"idx":[["BC","XB"]]},"deprecated":false,"name":"bulk_teloprimeV2"}}'
    )

    with tempfile.NamedTemporaryFile(delete=False) as f:
        f.write(json_string)

    lb_model = bam_utils.load_model(f.name)
    stored_model = model.LibraryModel.from_json_obj(json.loads(json_string))

    _compare_models(lb_model, stored_model)


def test_reverse_complement():
    d = {
        "GTTCTAGCGCTAGTATG": "CATACTAGCGCTAGAAC",
        "CTGCATAAAT": "ATTTATGCAG",
        "TTTCGCGCATATAG": "CTATATGCGCGAAA",
        "CCAA": "TTGG",
        "TTATATTTAT": "ATAAATATAA",
        "TANG": "CNTA",
        "SGYACGR": "YCGTRCS",
        "GCWATMCKVAABHGCD": "HGCDVTTBMGKATWGC",
    }

    for k, v_exp in d.items():
        v_obs = bam_utils.reverse_complement(k)
        v_obs_lc = bam_utils.reverse_complement(k.lower())

        assert v_obs == v_exp
        assert v_obs_lc == v_exp.lower()


@pytest.mark.slow
def test_generate_read_name_produces_no_collisions():
    movie_name = "m64020e_210000_000000"

    read_names = set()
    for zmw in range(8000000):
        for split_read_index in range(20):
            new_read_name = bam_utils.generate_read_name(
                movie_name, zmw, split_read_index
            )

            assert new_read_name not in read_names

            read_names.add(new_read_name)
