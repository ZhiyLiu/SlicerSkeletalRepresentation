# assume the folder is arranged by case id
# in each case folder, there are image file, ending with vtk, srep files(xml, vtps)
# output folder path
# this script will create subfolders in output folder for each case
import vtk, qt, ctk, slicer
import subprocess
import os
import re
import sys
import argparse

parser = argparse.ArgumentParser(description="Batch Structure Set Conversion")
parser.add_argument("-i", "--input-folder", dest="input_folder", metavar="PATH",
                    default="-", required=True,
                    help="Folder storing s-reps and surface meshes")
parser.add_argument("-s", "--initial-region", dest="initial_region_size", type=float,
                    default=0.01, required=False,
                    help="initial region size for NEWUOA")
parser.add_argument("-e", "--final-region", dest="final_region_size", type=float,
                    default=0.001, required=False, 
                    help="Use anatomy image as reference when converting structure set to labelmap")
parser.add_argument("-n", "--max-iter", dest="max_iter", type=int,
                    default=2000, required=False, 
                    help="Max iteration")
parser.add_argument("-wi", "--weight-image-match", dest="w_img_match", type=float,
                    default=0.004, required=False, 
                    help="Weight image match")
parser.add_argument("-wn", "--weight-normal-dev", dest="w_normal", type=float,
                    default=20.0, required=False,
                    help="Weight normal deviation")
parser.add_argument("-wg", "--weight-geo", dest="w_geo", type=float,
                    default=50, required=False, 
                    help="Weight geometric illegality")
parser.add_argument("-l", "--interp-level", dest="interp_level", type=int,
                        default=3, required=False,
                        help="Interpolation level")

parser.add_argument("-o", "--output-folder", dest="output_folder", metavar="PATH",
                    default=".", required=True,
                    help="Folder for output refined s-reps")

args = parser.parse_args(sys.argv[1:])

parent_dir = args.input_folder
output_dir = args.output_folder
init_size = args.initial_region_size
end_region = args.final_region_size
max_iter = args.max_iter
w_img = args.w_img_match
w_normal = args.w_normal
w_geo = args.w_geo
w_img = 0.00001
w_normal = 0.00001
w_geo = 10000
interp_level = args.interp_level
for patient_folder in os.listdir(parent_dir):
    case_folder = os.path.join(parent_dir, patient_folder)

    srep_file_path = os.path.join(case_folder, 'header.xml')
    img_file_path = None
    for filename in os.listdir(case_folder):
        if re.match(r"(.*)\.vtk", filename) != None:
            img_file_path = os.path.join(case_folder, filename)
            output_path = os.path.join(output_dir, patient_folder)
            if not os.path.exists(output_path):
                os.mkdir(output_path)
            refined_header_file = os.path.join(output_path, 'header.xml')
            if os.path.exists(refined_header_file):
                print(refined_header_file)
#            slicer.modules.skeletalrepresentationrefiner.logic().CLIRefine(srep_file_path, img_file_path, output_path, init_size, end_region, max_iter, w_img, w_normal, w_geo, interp_level)


print('Done!')