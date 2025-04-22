# **searchingforstrangeness**

## Building the Project

1. **Set up the environment:**

   ```bash
   apptainer shell \
          -B /cvmfs \
          -B /exp/uboone \
          -B /pnfs/uboone \
          -B /run/user \
          -B /etc/hosts \
          -B /etc/localtime \
          -s /bin/bash \
          --env UPS_OVERRIDE='-H Linux64bit+3.10-2.17' \
          /cvmfs/uboone.opensciencegrid.org/containers/uboone-devel-sl7
   source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh 
   setup ubooncode v08_00_00_80 -q e17:prof
   unsetup mrb
   setup mrb -o
   ```

   - This configures the environment using the `ubooncode` setup script with version `v08_00_00_89` and profile `e17:prof`, then sets up the `mrb` tool.

2. **Create a development area:**

   ```bash
   mkdir my_dev_area
   cd my_dev_area
   mrb newDev
   ```

   - Creates a new directory `my_dev_area` and initialises it as a development environment with `mrb`.

3. **Fetch the main repository:**

   ```bash
   mrb g ubana
   ```

   - Clones the `ubana` repository into the `srcs/ubana` directory using the `mrb` tool.


   cd srcs/ubana
   vim CMakeLists.txt 

   And modify the compiler flags so that it is this: 
   
cet_set_compiler_flags(DIAGS CAUTIOUS
  #WERROR
  NO_UNDEFINED
  ALLOW_DEPRECATIONS
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-expansion-to-defined -Wno-variadic-macros -Wno-pedantic 
  #-Wno-error=unused-variable 
)

i.e. comment out #WERROR and add -Wno-variadic-macros -Wno-pedantic 

4. **Clone additional repository:**

   ```bash
   cd srcs/ubana/ubana
   git clone https://github.com/njlane314/searchingforstrangeness.git
   ```

5. **Modify CMakeLists.txt:**

   ```bash
   # Add the following line to CMakeLists.txt to include searchingforstrangeness in the build
   echo "add_subdirectory(searchingforstrangeness)" >> CMakeLists.txt
   ```

   - Appends the `add_subdirectory(searchingforstrangeness)` line to `CMakeLists.txt` in `srcs/ubana/ubana`, integrating the cloned repository into the build process.

6. **Checkout the desired tag:**

   ```bash
   cd ../..  # Go back to srcs/ubana
   git checkout tags/v08.00.00
   ```

   - Returns to the `srcs/ubana` directory and checks out the specified tag `v08.00.00` for the `ubana` repository.

7. **Set up the build environment:**

   ```bash
   cd $MRB_TOP
   mrbsetenv
   ```

   - Moves to the top-level directory (defined by `$MRB_TOP`) and sets up the build environment with `mrbsetenv`.

8. **Build the project:**

   ```bash
   mrb i -j4
   ```

   - Builds the project using `mrb` with 4 parallel jobs (`-j4`).

## Processing Files

### **Local Processing**

To process ROOT files locally, you can manually run a series of commands to fetch, process, combine, and clean up files. Below are the key commands and an example of how to use them.

#### Key Commands

1. **Fetch a List of Files from a SAM Dataset**:
   ```bash
   files=$(samweb list-files defname:<sam_definition> | head -n <num_files>)
   ```
   - Retrieves a specified number of files from a SAM dataset.
   - **Example**: `files=$(samweb list-files defname:my_dataset | head -n 5)` fetches the first 5 files.

2. **Locate the File Path**:
   ```bash
   filedir=$(samweb locate-file <file> | grep -o '/pnfs/.*' | head -n 1)
   ```
   - Finds the physical directory path of a file.
   - **Example**: `filedir=$(samweb locate-file somefile.root | grep -o '/pnfs/.*' | head -n 1)`.

3. **Process a File with `lar`**:
   ```bash
   lar -c <fhicl_file> -s <filepath> -T <outputfile>
   ```
   - Runs the `lar` analysis tool on a file using a FHiCL configuration.
   - **Example**: `lar -c analysis.fcl -s /pnfs/uboone/some/path/somefile.root -T output_1.root`.

4. **Combine Output Files**:
   ```bash
   hadd -f <combined_output> <outputfiles>
   ```
   - Merges multiple ROOT files into one.
   - **Example**: `hadd -f combined_output.root output_*.root`.

5. **Clean Up Temporary Files**:
   ```bash
   rm <outputfiles>
   ```
   - Deletes temporary output files.
   - **Example**: `rm output_*.root`.

#### Example: Process 3 Files Locally
```bash
# Fetch 3 files from the dataset
files=$(samweb list-files defname:my_dataset | head -n 3)

# Process each file
counter=1
for file in $files; do
    filedir=$(samweb locate-file $file | grep -o '/pnfs/.*' | head -n 1)
    filepath="${filedir}/${file}"
    lar -c analysis.fcl -s $filepath -T output_$counter.root
    counter=$((counter + 1))
done

# Combine outputs
hadd -f combined_output.root output_*.root

# Clean up
rm output_*.root
```

## Configuration System

The framework uses FHiCL (Fermilab Hierarchical Configuration Language) for configuration.

### FHiCL Structure

The main configuration files include:

- **run_signal_selectionfilter.fcl**: Top-level configuration
  ```fcl
  process_name: SelectionSignalFilterProcess
  services: {
      TFileService: { fileName: "output.root" }
      # Detector configuration
      Geometry: @local::microboone_geo
      DetectorPropertiesService: @local::microboone_detproperties
      # Space charge correction
      SpaceCharge.EnableCorrSCE: true
  }
  physics: {
      filters: {
          selectionfilter: {
              module_type: SelectionFilter
              SelectionTool: @local::TruthSignalSelection
              EventType: "signal"  # Look for signal events only
              EventClassifier: @local::SharedEventClassifier
              AnalysisTools: {
                  neutrino: @local::NeutrinoAnalysis
                  pattern: @local::PatternAnalysis
                  wireimage: @local::WireImageAnalysis
                  # Other analysis tools...
              }
          }
      }
      trigger_paths: [ e1 ]
      e1: [ selectionfilter ]
  }
  ```

- **selectionconfig.fcl**: Defines tool configurations
  ```fcl
  # Event classifier configuration
  SharedEventClassifier: {
      SignatureTools: {
          leptonic: @local::MuonSignature
          hadronic: @local::LambdaSignature  # Change this to search for different particles
      }
  }
  
  # Selection tool definitions
  TruthSignalSelection: {
      tool_type: "TruthSignalSelection"
      EventClassifier: @local::SharedEventClassifier
      EventType: "signal"
  }
  ```

1. **Event Selection**:
   - Modify `EventType` to select "signal" or "background" events, avoiding double counting signal events that can appear in background samples. 
   - Configure `SignatureTools` section to define what combinations of particles constitute your chosen signal definition. 

2. **Analysis Tools**:
   - Configure tools from the `AnalysisTools` section.
   - Each tool produces a different subset of the output data, and can be removed or added. 


### **Grid Submission**

To submit jobs to the grid, you can manually run commands to package your code, authenticate, and submit jobs using `project.py`. Below are the key commands and an example.

The project uses XML files to define grid jobs, modify these where necessary:

```xml
<job>
<project name="&name;">
  <numevents>-1</numevents>
  <resource>DEDICATED,OPPORTUNISTIC,OFFSITE</resource>
  <larsoft>
    <tag>&release;</tag>
    <qual>e17:prof</qual>
    <local>/pnfs/uboone/resilient/users/nlane/NeutralKaon/tarballs/StrangenessCode.tar</local>
  </larsoft>
  <stage name="analyse">
    <inputdef>prod_strange_resample_fhc_run2_fhc_reco2_reco2</inputdef>
    <fcl>run_signal_selectionfilter.fcl</fcl>
    <outdir>/pnfs/uboone/scratch/users/nlane/kaon_dl/&release;/&name;/out</outdir>
    <memory>4000</memory>
    <disk>20GB</disk>
    <jobsub>--expected-lifetime=24h</jobsub>
  </stage>
</project>
</job>
```

#### Key Commands

1. **Create a Tarball**:
   ```bash
   make_tar_uboone.sh <tarball_name>
   ```
   - Packages your code into a tarball.
   - **Example**: `make_tar_uboone.sh my_tarball.tar`.

2. **Copy the Tarball**:
   ```bash
   cp -f <tarball_name> <tarball_dest>
   ```
   - Copies the tarball to a resilient directory.
   - **Example**: `cp -f my_tarball.tar /pnfs/uboone/resilient/users/myuser/`.

3. **Authenticate for Grid Access**:
   ```bash
   htgettoken -a <vault_server> -i <experiment>
   ```
   - Obtains an authentication token.
   - **Example**: `htgettoken -a htvaultprod.fnal.gov -i uboone`.

4. **Submit Jobs to the Grid**:
   ```bash
   project.py --xml <xml_config_file> --stage <stage> --submit
   ```
   - Submits jobs using a configuration file.
   - **Example**: `project.py --xml config.xml --stage analyse --submit`.

5. **Clean and Retry**:
   ```bash
   project.py --xml <xml_config_file> --stage <stage> --clean
   project.py --xml <xml_config_file> --stage <stage> --submit
   ```
   - Cleans up and retries submission.
   - **Example**: `project.py --xml config.xml --stage analyse --clean` followed by `project.py --xml config.xml --stage analyse --submit`.

#### Example: Submit Jobs to the Grid
```bash
# Create and copy tarball 
make_tar_uboone.sh my_tarball.tar
cp -f my_tarball.tar /pnfs/uboone/resilient/users/myuser/

# Authenticate
htgettoken -a htvaultprod.fnal.gov -i uboone

# Submit jobs
project.py --xml config.xml --stage analyse --submit

# If it fails, clean and retry
project.py --xml config.xml --stage analyse --clean
project.py --xml config.xml --stage analyse --submit
```

## General Commands

### **Environment and Container Management**

- **`apptainer shell [options] /path/to/container`**  
    Launches an interactive shell within an Apptainer container, providing an isolated environment for running applications.  
    - **Example**: 
        ```bash
        apptainer shell \
                -B /cvmfs \
                -B /exp/uboone \
                -B /pnfs/uboone \
                -B /run/user \
                -B /etc/hosts \
                -B /etc/localtime \
                -s /bin/bash \
                --env UPS_OVERRIDE='-H Linux64bit+3.10-2.17' \
                /cvmfs/uboone.opensciencegrid.org/containers/uboone-devel-sl7
        ```

- **`source setup.sh`**  
    Executes a script to configure the shell environment, often setting variables or loading dependencies.  
    - **Example**: `source /setup.sh` loads necessary configurations.
    - **Script**: 
        ```bash
        source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
        setup uboonecode v08_00_00_82 -q "e17:prof"
        ```

- **`setup <product> <version> -q <qualifiers>`**  
    Configures the environment for a specific software product and version using UPS (Unix Product Support).  
    - **Example**: see above

---

### **Authentication and Grid Access**

- **`kx509`**  
  Converts a Kerberos ticket into an X.509 certificate for secure authentication to access grid resources.  
  - **Example**: Run `kx509` after `kinit` to obtain a certificate.

- **`voms-proxy-init -noregen -voms <virtual_organization>`**  
  Generates a proxy certificate for grid authentication tied to a virtual organization.  
  - **Example**: `voms-proxy-init -noregen -voms fermilab:/fermilab/uboone/Role=Analysis` creates a proxy for Fermilab.

---

### **File and Data Management**

- **`find /path/to/directory -type f -name "*.root"`**  
  Searches for files matching a pattern (e.g., `*.root` files).  

- **`hadd -f output.root file1.root file2.root ...`**  
  Merges multiple ROOT files into one.  
  - **Example**: `hadd -f combined.root run1.root run2.root` merges two files.

- **`tar -cvf tarball_name.tar /path/to/code`**  
  Creates a tarball of a directory or files.  

- **`samweb list-files defname:<sam_definition>`**  
  Lists files in a SAM dataset.  

- **`samweb locate-file <file>`**  
  Finds the physical location of a SAM-managed file.  

- **`samweb count-files defname:<definition_name>`**  
  Counts files in a SAM dataset.  

- **`samweb prestage-dataset --defname=<definition_name>`**  
  Stages a SAM dataset for fast access.  

- **`samweb create-definition my_subset "defname:make_lambda_overlay_nohadrons_reco2_reco2 with limit 100"`**
   Creates limited dataset 

---

### **Job Submission and Monitoring**

- **`jobsub_q $USER`**  
  Queries the status of a user's grid jobs.  

- **`project.py --xml <xml_file> --stage <stage> --submit`**  
  Submits jobs to the grid, automates workflow with `*.xml` configurations.  

---

### **Data Processing and Analysis**

- **`lar -c <fhicl_file> -s <input_file> -T <output_file>`**  
  Runs LArSoft with a FHiCL configuration.  

- **`lar -c eventdump.fcl -s <your.root> | grep recob::Wire`**  
  Dumps event information from a `*.root` file and searches for `recob::Wire` objects.

---

### **Machine Learning and Library Setup**

- **`setup libtorch <version> -q <qualifiers>`**  
  Configures libtorch for machine learning tasks.

### **Grid Submission**

## NuMI Input Samdefs

- **New_NuMI_Flux_Run_1_FHC_Pandora_Reco2_reco2_reco2      # Beam background**
- **prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2    # EXT/Beam-off background**
- **prod_strange_resample_fhc_run2_fhc_reco2_reco2         # Enriched strangeness**
- **make_lambda_overlay_nohadrons_reco2_reco2              # Lambda no-hadrons**
- **prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_v2_sample0              # Background**
- **cthorpe_make_k0s_events_numi_rhc_reco2_REAL_reco2_reco2                  # Kaons**
- **cthorpe_prod_extnumi_mcc9_v08_00_00_45_run3_run3b_reco2_all_reco2_pt1    # EXT**
- **cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2           # Hyperons**




source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
setup larsoft v10_04_07 -q e26:prof
setup libtorch v2_1_1a -q e26
setup mrb

cd <ProductionDir>
mrb newDev
source localProducts_*/setup

cd srcs
mrb g ubana
cd ubana
git checkout tags/v08_00_00_82

cd $MRB_TOP
mrbsetenv
mrb i -j4






cd [nlane@uboonegpvm03 ~]$ cd /exp/uboone/app/users/nlane/production/
[nlane@uboonegpvm03 production]$ 
[nlane@uboonegpvm03 production]$ 
[nlane@uboonegpvm03 production]$ 
[nlane@uboonegpvm03 production]$ 
[nlane@uboonegpvm03 production]$ apptainer shell           -B /cvmfs           -B /exp/uboone           -B /pnfs/uboone           -B /run/user           -B /etc/hosts           -B /etc/localtime           -s /bin/bash           --env UPS_OVERRIDE='-H Linux64bit+3.10-2.17'           /cvmfs/uboone.opensciencegrid.org/containers/uboone-devel-sl7
Apptainer> 
Apptainer> 
Apptainer> 
Apptainer> 
Apptainer> ls
build_instructions.txt	deprecated  strangeness_mcc10  strangeness_mcc9
Apptainer> rm -rf strangeness_mcc9 
Apptainer> mkdir strangeness_mcc9 
Apptainer> 
Apptainer> 
Apptainer> 
Apptainer> cd strangeness_mcc
bash: cd: strangeness_mcc: No such file or directory
Apptainer> 
Apptainer> 
Apptainer> cd strangeness_mcc9 
Apptainer> 
Apptainer> 
Apptainer> 
Apptainer> 
Apptainer> mrb newDev -q e17:prof
bash: mrb: command not found
Apptainer> setup mrb
bash: setup: command not found
Apptainer> source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh 
Setting LC_ALL=C
Scientific Linux 7.9 (Nitrogen)
Setting up larsoft UPS area... /cvmfs/larsoft.opensciencegrid.org/products
Setting up uboone UPS area... /cvmfs/uboone.opensciencegrid.org/products
Setting up fermilab common UPS area... /cvmfs/fermilab.opensciencegrid.org/products/common/db
Apptainer> mrb newDev -q e17:prof

building development area for larsoft  -q e17:prof


The following configuration is defined:
  The top level directory is .
  The source code directory will be under .
  The build directory will be under .
  The local product directory will be under .

MRB_BUILDDIR is /exp/uboone/app/users/nlane/production/strangeness_mcc9/build_slf7.x86_64
MRB_SOURCE is /exp/uboone/app/users/nlane/production/strangeness_mcc9/srcs 
INFO: cannot find larsoft//releaseDB/base_dependency_database
      or larsoftcode//releaseDB/base_dependency_database
      mrb checkDeps and pullDeps will not have complete information

IMPORTANT: You must type
    source /exp/uboone/app/users/nlane/production/strangeness_mcc9/localProducts_larsoft__e17_prof/setup
NOW and whenever you log in

Apptainer> source /exp/uboone/app/users/nlane/production/strangeness_mcc9/localProducts_larsoft__e17_prof/setup

MRB_PROJECT=larsoft
MRB_PROJECT_VERSION=
MRB_QUALS=e17:prof
MRB_TOP=/exp/uboone/app/users/nlane/production/strangeness_mcc9
MRB_SOURCE=/exp/uboone/app/users/nlane/production/strangeness_mcc9/srcs
MRB_BUILDDIR=/exp/uboone/app/users/nlane/production/strangeness_mcc9/build_slf7.x86_64
MRB_INSTALL=/exp/uboone/app/users/nlane/production/strangeness_mcc9/localProducts_larsoft__e17_prof

PRODUCTS=/exp/uboone/app/users/nlane/production/strangeness_mcc9/localProducts_larsoft__e17_prof:/cvmfs/uboone.opensciencegrid.org/products:/cvmfs/larsoft.opensciencegrid.org/products:/cvmfs/fermilab.opensciencegrid.org/products/common/db

Apptainer> pwd
/exp/uboone/app/users/nlane/production/strangeness_mcc9
Apptainer> mrb g ubana
Cloning into 'ubana'...
remote: Enumerating objects: 66440, done.
remote: Counting objects: 100% (88/88), done.
remote: Compressing objects: 100% (73/73), done.
remote: Total 66440 (delta 39), reused 36 (delta 15), pack-reused 66352 (from 1)
Receiving objects: 100% (66440/66440), 27.73 MiB | 39.11 MiB/s, done.
Resolving deltas: 100% (49193/49193), done.
NOTICE: Adding ubana to CMakeLists.txt file
Apptainer> pwd
/exp/uboone/app/users/nlane/production/strangeness_mcc9
Apptainer> ls
build_slf7.x86_64  localProducts_larsoft__e17_prof  srcs
Apptainer> cd srcs/
Apptainer> 
Apptainer> 
Apptainer> ls
CMakeLists.txt  ubana
Apptainer> cd ubana/          
Apptainer> git checkout tags/v08_00_00_82
Note: switching to 'tags/v08_00_00_82'.

You are in 'detached HEAD' state. You can look around, make experimental
changes and commit them, and you can discard any commits you make in this
state without impacting any branches by switching back to a branch.

If you want to create a new branch to retain commits you create, you may
do so (now or later) by using -c with the switch command. Example:

  git switch -c <new-branch-name>

Or undo this operation with:

  git switch -

Turn off this advice by setting config variable advice.detachedHead to false

HEAD is now at e0379f989 v08_00_00_82
Apptainer> cd ../..
Apptainer> ls
build_slf7.x86_64  localProducts_larsoft__e17_prof  srcs
Apptainer> pwd
/exp/uboone/app/users/nlane/production/strangeness_mcc9
Apptainer> mrbsetenv
The working build directory is /exp/uboone/app/users/nlane/production/strangeness_mcc9/build_slf7.x86_64
The source code directory is /exp/uboone/app/users/nlane/production/strangeness_mcc9/srcs
----------- check this block for errors -----------------------
----------------------------------------------------------------
Apptainer> mrb b