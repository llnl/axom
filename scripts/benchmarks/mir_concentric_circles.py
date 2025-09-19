import os, sys
import argparse
import math

def flux_run(np):
  return f"flux run -n {np}"

def srun(np):
  return f"srun -n {np}"

def lrun(np):
  return f"lrun -n {np}"

# Parallel options use 4 ranks per node since that is how many GPUs there are.
runs = {
  "build-rzansel-blueos_3_ppc64le_ib_p9-clang@10.0.1.2_cuda-release" : {"policies":["seq", "omp", "cuda"], "launch":lrun},
  "build-rzwhippet-toss_4_x86_64_ib-clang@14.0.6-release" : {"policies":["seq", "omp"], "launch":srun},
  "build-rzwhippet-toss_4_x86_64_ib-gcc@10.3.1-release" : {"policies":["seq", "omp"], "launch":srun},
  "build-rzwhippet-toss_4_x86_64_ib-intel@2022.1.0-release" : {"policies":["seq", "omp"], "launch":srun},
  "build-rzvernal-toss_4_x86_64_ib_cray-rocmcc@6.2.1_hip-release" : {"policies":["seq", "hip"], "launch":srun},
  "build-rzvernal-toss_4_x86_64_ib_cray-clang@17.0.0_hip-release" : {"policies":["seq", "hip"], "launch":srun},
  "build-rzadams-toss_4_x86_64_ib_cray-cce@18.0.0_hip-release" :{"policies":["seq", "hip"], "launch":flux_run},
  "build-rzadams-toss_4_x86_64_ib_cray-rocmcc@6.2.1_hip-release" :{"policies":["seq", "hip"], "launch":flux_run}
}

def generate(params):
  """
  Generate scripts to run the program to create timings.
  """
  method = params["method"]

  for r in runs:
    if not os.path.exists(r):
      print(f"Skipping {r}")
      continue
    filename = os.path.join(r, "run_concentric_circles.bash")

    f = open(filename, "wt")
    f.write("#!/bin/bash\n\n")
    f.write("CONCENTRIC_CIRCLES=./examples/mir_concentric_circles\n")
    f.write("CONCENTRIC_CIRCLES_MPI=./examples/mir_concentric_circles_mpi\n\n")

    dimension = params["dimension"]
    trials = params["trials"]
    for s in params["sizes"]:
      f.write(f"# Size {s}\n")
      if len(params["parallel"]) > 0:
        # parallel
        for np in params["parallel"]:
          launch = runs[r]["launch"](np)
          for policy in runs[r]["policies"]:
            f.write(f'echo "Running {launch} $CONCENTRIC_CIRCLES_MPI --gridsize {s} --numcircles 5 --policy {policy} --method {method} --dimension {dimension} --trials {trials} --disable-write"\n')
            f.write(f'{launch} $CONCENTRIC_CIRCLES_MPI --gridsize {s} --numcircles 5 --policy {policy} --method {method} --dimension {dimension} --trials {trials} --disable-write > result_{policy}_np{np}_s{s}.txt\n\n')
      else:
        # serial
        for policy in runs[r]["policies"]:
          f.write(f'echo "Running --gridsize {s} --numcircles 5 --policy {policy} --method {method} --dimension {dimension} --trials {trials} --disable-write"\n')
          f.write(f'$CONCENTRIC_CIRCLES --gridsize {s} --numcircles 5 --policy {policy} --method {method} --dimension {dimension} --trials {trials} --disable-write > result_{policy}_s{s}.txt\n\n')

    f.close()
    os.chmod(filename, 0o700)
    print(f"Wrote {filename}")

def read_timings(filename, searchKey):
  """
  Read a single Caliper log and pull out the timings we want.
  """
  retval = "" # no data
  try:
    lines = open(filename, "rt").readlines()
    print(f"Reading {filename}")
    for line in lines:
      pos = line.find(searchKey)
      if pos != -1:
        toks = [x for x in line.split() if x != '']
        #print(toks)
        try:
          retval = float(toks[2]) # timings (I)
        except ValueError:
          # The line had our searchKey but was not a Caliper line.
          continue
        break
  except:
    print(f"Error reading {filename}")
    pass
  return retval

def make_columns(params):
  """
  Read the available timing files and make columns of data.
  """
  # I was measuring runMIR before, which includes device transfers and HDF5 output.
  # I don't want that stuff in the timings because it is not relevant to timing MIR
  # and it can vary wildly.

  def is_selected(selections, name):
    """
    Checks whether a name is being selected.
    """
    selected = True
    if len(selections) > 0:
      selected = False
      for s in selections:
        if name.find(s) != -1:
          selected = True
          break
    return selected 

  def average(value, trials):
    try:
      avg = str(float(value) / float(trials))
    except ValueError:
      avg = value
    return avg

  # Measure just the MIR algorithm
  searchKey = "EquizAlgorithm"
  if params["method"] == "elvira":
    searchKey = "ElviraAlgorithm"

  columns = []
  # Add NumZones column (either square or cube of s, depending on dimension)
  sc = ["NumZones"]
  if params["dimension"] == 3:
    for s in params["sizes"]:
      sc.append(s*s*s)
  else:
    for s in params["sizes"]:
      sc.append(s*s)
  columns.append(sc)

  trials = params["trials"]

  # Gather data.
  for r in sorted(runs.keys()):
    # Make sure r exists
    if not os.path.exists(r):
      continue

    # Gather individual data files
    for policy in runs[r]["policies"]:
      buildname = r[6:-8]
      if len(params["parallel"]) > 0:
        for np in params["parallel"]:
          name = f"{buildname} {policy.upper()} NP={np}"
          # Check if this name is selected
          if not is_selected(params["selections"], name):
             continue
          data = [name]
          for s in params["sizes"]:
            filename = os.path.join(r, f"result_{policy}_np{np}_s{s}.txt")
            value = read_timings(filename, searchKey)
            # Timings contain #trials values so average the value.
            avgValue = average(value, trials)
            data.append(avgValue)
          columns.append(data)
      else:
        name = f"{buildname} {policy.upper()}"
        data = [name]
        # Check if this name is selected
        if not is_selected(params["selections"], name):
           continue
        for s in params["sizes"]:
          filename = os.path.join(r, f"result_{policy}_s{s}.txt")
          value = read_timings(filename, searchKey)
          # Timings contain #trials values so average the value.
          avgValue = average(value, trials)
          data.append(avgValue)
        columns.append(data)
  return columns

def make_csv(params, outputfile):
  """
  Read the available timing files and assemble a CSV file.
  """
  columns = make_columns(params)

  # Write data
  f = open(outputfile, "wt")
  nrows = len(columns[0])
  for i in range(nrows):
    rowdata = [str(c[i]) for c in columns]
    line = ",".join(rowdata)
    f.write(f"{line}\n")
  f.close()

def seriesName(name):
  policies = {"SEQ" : "Serial", "OMP": "OpenMP", "CUDA" : "CUDA", "HIP" : "HIP"}
  newName = name[:name.find("-")]
  for p in policies.keys():
    if name[-len(p) - 1:] == " " + p:
      newName = newName + " " + policies[p]
      break
  return newName

def lineProps(name):
  """
  Make line properties when given a series name.
  """
  hostColor = {"rzansel" : "g", "rzwhippet" : "r", "rzvernal" : "o", "rzadams" : "b"}
  policyStyle = {"seq" : "--", "omp": ":", "cuda" : "-", "hip" : "-"}
  policyMark = {"seq" : "o", "omp": "s", "cuda" : "^", "hip" : "^"}
  color = "b"
  style = "-"
  mark = "o"
  for h in hostColor.keys():
     if name.find(h) != -1:
       color = hostColor[h]
       for p in policyStyle.keys():
         if name[-len(p) - 1:].lower() == " " + p:
           style = policyStyle[p]
           mark = policyMark[p]
           break
       break
  return color, style, mark

def plot(params):
  """
  Read the available timing files and plot them.
  """
  def make_series(col1, col2):
    n = len(col1)
    x = []
    y = []
    for i in range(n):
      try:
        xv = float(col1[i])
        yv = float(col2[i])
        x.append(xv)
        y.append(yv)
      except ValueError:
        # We probably did not get data. Stop processing.
        break
    return x,y

  columns = make_columns(params)

  import matplotlib.pyplot as plt
  for c in range(1, len(columns)):
    x, y = make_series(columns[0][1:], columns[c][1:])
    if len(x) > 0:
      color, style, mark = lineProps(columns[c][0])
      plt.plot(x, y, marker=mark, linestyle=style, color=color, linewidth=2., label=seriesName(columns[c][0]))

  dimension = params["dimension"]
  method = params["method"]

  # Add labels and title
  plt.xlabel('Number of Zones', fontsize=24)
  plt.ylabel('Time (s)', fontsize=24)
  plt.title(f'{dimension}D MIR Timings ({method})', fontsize=28)
  xlabels = make_series(columns[0][1:], columns[0][1:])[0]
  plt.xticks(ticks=xlabels, labels=xlabels, fontsize=18)
  plt.yticks(fontsize=18)

  # Set x-axis to logarithmic scale
  plt.xscale('log')

  # Set y-axis to logarithmic scale
  plt.yscale('log')

  # Add a legend
  plt.legend(fontsize=18)

  try:
    plt.tight_layout()
  except ValueError:
    print(f"There was an error, probably because the benchmark is still running. {columns}")
    raise

  # Show the plot
  plt.grid(True)
  plt.show()

def get_params():
  """
  Get arguments from the command line and return a dictionary.
  """
  parser = argparse.ArgumentParser(description="Parse args")
  parser.add_argument(
    "--parallel",
    type=str,
    help="Comma-separated list of integer values",
    required=False
  )

  parser.add_argument(
    "--generate",
    action="store_true",
    help="Boolean flag to indicate whether to generate scripts",
    required=False
  )

  parser.add_argument(
    "--plot",
    action="store_true",
    help="Boolean flag to indicate whether to plot",
    required=False
  )

  parser.add_argument(
    "--method",
    type=str,
    help="MIR method to use (e.g., 'equiz', 'elvira')",
    required=False
    )

  parser.add_argument(
    "--dimension",
    type=int,
    help="Mesh dimension to generate (2 or 3)",
    required=False
    )

  parser.add_argument(
    "--select",
    type=str,
    help="Comma-separated list of strings to select specific plot data",
    required=False
  )

  parser.add_argument(
    "--trials",
    type=int,
    help="Number of times to run MIR",
    required=False
    )

  args = parser.parse_args()

  # Convert the comma-separated string into a tuple of integers
  parallel = ()
  if args.parallel is not None:
    try:
      parallel = tuple(int(value) for value in args.parallel.split(","))
    except ValueError:
      raise argparse.ArgumentTypeError("All values in --parallel must be integers.")

  selections = ()
  if args.select is not None:
    selections = tuple(args.select.split(","))

  params = {}
  params["parallel"] = parallel
  if args.method is not None:
    params["method"] = args.method
  else:
    params["method"] = "elvira"

  if args.generate:
    params["generate"] = args.generate
  else:
    params["generate"] = False

  if args.plot:
    params["plot"] = args.plot
    params["generate"] = False
  else:
    params["plot"] = False

  params["selections"] = selections
  if args.dimension is not None:
    params["dimension"] = args.dimension
  else:
    params["dimension"] = 2

  if args.trials is not None:
    params["trials"] = max(1, args.trials)
  else:
    params["trials"] = 1

  # Generate some sizes.
  sides = (50, 100, 200, 500, 1000, 1500, 2000, 4000, 8000)
  if params["dimension"] == 3:
    # Figure out the 3D side size we need for about the same zone counts as 2D.
    s = []
    for side in sides:
      s.append(int(math.ceil(math.pow(side * side, 1./3.))))
    params["sizes"] = tuple(s)
  else:
    params["sizes"] = sides

  return params

def main():
  params = get_params()

  if params["generate"]:
    print("Generating files...")
    generate(params)
  elif params["plot"]:
    print("Plotting files...")
    plot(params)
  else:
    print("Making CSV...")
    make_csv(params, "concentric_circle_timings.csv")

if __name__ == "__main__":
  main()
