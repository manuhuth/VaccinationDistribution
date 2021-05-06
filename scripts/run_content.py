#!/usr/bin/env python
"""This script compiles the paper."""
import argparse
import shutil
import os

ROOT = os.environ["PROJECT_ROOT"]


def compile_latex_document(dirname=None):
    cwd = os.getcwd()
    if dirname:
        os.chdir(dirname)

        literature = "bibtex"
        if "slides" in dirname:
            literature = "biber"

    [
        os.system(type_ + " main")
        for type_ in ["pdflatex", literature, "pdflatex", "pdflatex"]
    ]

    os.chdir(cwd)


if __name__ == "__main__":

    parser = argparse.ArgumentParser("Create content")

    parser.add_argument("-p", "--paper", action="store_true", help="compile paper")

#    parser.add_argument("-s", "--slides", action="store_true", help="compile slides")

    args = parser.parse_args()

    os.chdir(ROOT)

    if args.paper:
        compile_latex_document(ROOT + "/paper")
        shutil.copy("paper/main.pdf", "main.pdf")

#    if args.slides:
#        compile_latex_document(ROOT + "")
#        shutil.copy("", "slides.pdf")
