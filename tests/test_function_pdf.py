import basicsynbio as bsb
import pytest

from .test_fixtures import small_build_example


def test_pdf_is_built_of_expected_size(small_build_example):
    import PyPDF2
    import os

    build = small_build_example
    pdffilename = bsb.pdf_instructions(build)
    file = open(pdffilename, "rb")
    fileReader = PyPDF2.PdfFileReader(file)
    assert fileReader.numPages == 4
    os.remove(pdffilename)


def test_correct_template_text_written_to_pdf(small_build_example):
    import PyPDF2
    import os

    build = small_build_example
    pdffilename = bsb.pdf_instructions(build)
    file = open(pdffilename, "rb")
    fileReader = PyPDF2.PdfFileReader(file)
    pageObj = fileReader.getPage(0)
    firstPageText = pageObj.extractText()
    assert "Materials" in firstPageText
    assert "Ambion AM10050 (Thermo)" in firstPageText
    os.remove(pdffilename)


def test_assembly_build_objects_are_written_to_pdf(small_build_example):
    import PyPDF2
    import os

    build = small_build_example
    pdffilename = bsb.pdf_instructions(build)
    file = open(pdffilename, "rb")
    fileReader = PyPDF2.PdfFileReader(file)
    allpdftext = ""
    for i in range(fileReader.numPages):
        pageObj = fileReader.getPage(i)
        pageText = pageObj.extractText()
        allpdftext += pageText
    processed_text = "".join(allpdftext.split())
    assert "BASIC_sfGFP_CDS" in processed_text
    assert "Second_Assembly_With_26" in processed_text
    os.remove(pdffilename)
