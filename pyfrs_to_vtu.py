import os
import subprocess
import argparse
import tempfile

ap = argparse.ArgumentParser()
ap.add_argument("mesh-file", metavar="mesh_f", type=str, help="Mesh file for creating .vtu output.")
ap.add_argument("-d", required=False, default=4, help="DIVISOR\nsets the level to which high order elements are divided; output is linear between nodes, so increased resolution may be required", type=int)
ap.add_argument("-g", "--gradients", help="compute gradients.", action='store_true')
ap.add_argument("-p", "--precision", required=False, default="single", help="output number precision; defaults to single", choices=['single', 'double'])
ap.add_argument("-i", "--in", required=False, type=str, help="Folder Containg .pyfrs files.", default="./")
ap.add_argument("-o", "--out", required=False, type=str, help="Folder to output files to.", default="./")
ap.add_argument("-t", required=False, type=str, help="Output file type; defaults to .vtu", default=".vtu")
ap.add_argument("pvd_path", type=str, help="Path to write .pvd file to.")

args = vars(ap.parse_args())

print(args)

processes = []
input_files  = os.listdir(args['in'])
input_files  = [i for i in input_files if i.endswith(".pyfrs")]

output_files = [i.replace(".pyfrs", args['t']) for i in input_files]

input_files  = [os.path.join(args['in'], i) for i in input_files]
output_files = [os.path.join(args['out'], i) for i in output_files]

base_command = ["pyfr","export"]

base_command.extend(["-d", str(args['d'])])
base_command.extend(["-p", str(args['precision'])])

if args['gradients']:
    base_command.extend(["-g"])

base_command.extend([args['mesh-file']])

for i_fle, o_fle in zip(input_files, output_files):
    full_command = base_command + [i_fle, o_fle]
    f = tempfile.TemporaryFile()
    p = subprocess.Popen(full_command, stdout=f)
    processes.append((p,f))

for p,f in processes:
    p.wait()
    f.seek(0)
    print("\n".join(f.readlines()))
    f.close()

pvd_path = args['pvd_path']

if not pvd_path.endswith(".pvd"):
    pvd_path += ".pvd"

pvd_start = """<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>"""

pvd_end = """  
  </Collection>
</VTKFile>
"""

pvd_entry_line = "<DataSet timestep=\"{tstep}\" file=\"{fle_path}\"/>"

timesteps = [float(i.split("-")[-1].replace(args['t'], "")) for i in output_files]
ordered_out = sorted(zip(timesteps, output_files), key=lambda x: x[0])

with open(os.path.join(args['out'], args['pvd_path']), 'w') as fle:
    fle.write(pvd_start)
    for tstep, fle_path in ordered_out:
        fle_path = fle_path.split("/")[-1]
        fle.write("\n  " + pvd_entry_line.format(tstep=tstep, fle_path=fle_path))
    fle.write(pvd_end)
