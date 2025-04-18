# **searchingforstrangeness**

_This [searchingforstrangeness](https://github.com/njlane314/searchingforstrangeness) code is a [Pandora](https://github.com/PandoraPFA/larpandora) neutrino selection framework, built with [LArSoft](https://github.com/LArSoft/larsoft), using the design implemented in [searchingfornues](https://github.com/ubneutrinos/searchingfornues)._

_The framework loads a candidate neutrino slice given by Pandora, runs a single selection tool for each event, and a series of analysis tools for each slice. The selection and analysis tools are configurable at run, and the output is to a single root file. This implementation is specifically designed to search for rare neutrino processes that give distinct topologies or identifiable pattern signatures within the MicroBooNE detector; these signatures are identified using a visual deep-learning network that processes each of the wire planes independently._

## Code Structure

### Core Components

- **SelectionFilter**: Main module that processes events and applies a selection, filtering events. 
   - `_selectionTool`: A tool that determines if an event passes a defined selection. 
   - `_analysisToolsVec`: A collection of analysis tools applied to each slice, and event. 
   - `_eventClassifier`: Categorises events by the configurable signal defintion and event type. 

- **EventClassifier**: Determines if an event is signal or background.
   - `_pattern`: Collection of event signatures defined by signature tools. 
   - `_signatureToolsVec`: Tools that define specific particle signatures from their decay topology.
   - `_clarityToolsVec`: Tools that evaluate the clarity of event topologies in the detector.

- **LabelClassifier**: Assigns particle-level truth classifications, basde on NuGraph definition. 
   - `_gamma_threshold`: Energy threshold for photon classification.
   - `_hadron_threshold`: Momentum threshold for hadron classification.

- **ImageProcessor**: Creates image representations of detector raw data.
   - `ImageProperties`: Defines dimensions and scaling of images.
   - `Image`: Contains pixel data and metadata for a single image.

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

   - This configures the environment using the `ubooncode` setup script with version `v08_00_00_80` and profile `e17:prof`, then sets up the `mrb` tool.

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
