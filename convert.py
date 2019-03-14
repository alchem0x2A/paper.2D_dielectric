import re
import os, os.path
import shutil
import glob
import subprocess
import multiprocessing
from multiprocessing import Pool
import warnings



# color definition
class TColors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def convert_pdf(infile, outdir="./img"):
    # Convert svg to pdf
    base_name = os.path.basename(infile)
    infile = os.path.join(PATH, base_name)
    outfile = os.path.join(os.path.abspath(outdir), base_name)
    outfile = outfile.replace(".svg", ".png")
    print(outfile)
    program = "/Applications/Inkscape.app/Contents/Resources/bin/inkscape"
    params = ["--without-gui", "--export-area-page",
    ]
    io = ["--file={}".format(infile),
          "--export-png={}".format(outfile),
          "--export-dpi=600"]
    success = subprocess.call([program, *params, *io])
    
    if success != 0:
        warnings.warn(TColors.FAIL + "File {} cannot be converted!".format(infile) + TColors.ENDC)
    else:
        print(TColors.OKGREEN +
              "File {} converted successfully on thread {}.".format(infile, multiprocessing.current_process())
              +TColors.ENDC)
    program = "convert"
    infile = outfile
    outfile = infile.replace(".png", ".pdf")
    io = [infile, outfile]
    success = subprocess.call([program, *io])
    if success != 0:
        warnings.warn(TColors.FAIL + "File {} cannot be converted!".format(infile) + TColors.ENDC)
    else:
        print(TColors.OKGREEN +
              "File {} converted to pdf on thread {}.".format(infile, multiprocessing.current_process())
              +TColors.ENDC)

# Convert all the pdf files
PATH = os.path.abspath("./img")

if __name__ == "__main__":
    import sys
    try:
        img_path = sys.argv[1]
    except IndexError:
        img_path = None
    file_list = []
    # TeX process
    # for ifile in glob.glob("*.tex"):
        # tex_process(ifile)
        # print(TColors.OKBLUE + "Converted TeX file: {}".format(ifile) + TColors.ENDC)

    # PDF Process 
    for ifile in glob.glob(os.path.join(PATH, "*.svg")):
        if img_path is None:
            file_list.append((ifile))
        else:
            file_list.append((ifile, img_path))
    N_cores = multiprocessing.cpu_count()
    # multicore
    with Pool(N_cores) as p:
        p.map(convert_pdf, file_list)
