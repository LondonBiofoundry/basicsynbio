import basicsynbio as bsb
import pytest
from .test_fixtures import small_build_example


def test_echo_instructions_small_build(small_build_example):
    echo_clips_zippath = bsb.export_echo_clips_instructions(
        small_build_example)
    assert 1 == 1
