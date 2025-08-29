# **searchingforstrangeness**

## Building the Project

Before building and running the inference utilities, ensure that required Python
packages such as `h5py` are available in the runtime environment.

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
   ```
    - Setup container and source dependencies. 

2. **Create a development area:**

   ```bash
   cd /exp/uboone/app/users/$USER/production/
   mkdir strangeness
   cd strangeness
   mrb newDev -q e17:prof
   source /exp/uboone/app/users/$USER/production/strangenes/localProducts_*/setup
   ```

   - Creates a new directory `strangeness` and initialise it as a development environment with `mrb`.

3. **Fetch the main repository:**

   ```bash
   mrb g ubana
   ```

   - Clones the `ubana` repository into the `srcs/ubana` directory using the `mrb` tool.

   ```bash
   cd srcs/ubana  # Go to srcs/ubana
   git checkout tags/v08_00_00_82
   ```
   
   ```bash
   cd srcs/ubana
   vim CMakeLists.txt
   # Modify the compiler flags as such 
   cet_set_compiler_flags(DIAGS CAUTIOUS
   #WERROR
   NO_UNDEFINED
   ALLOW_DEPRECATIONS
   EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-expansion-to-defined -Wno-variadic-macros -Wno-pedantic 
   #-Wno-error=unused-variable 
   )
   ```

   - Modify the compiler flags to disable treat warning as errors, this is needed for the LibTorch library. 

4. **Clone additional repository:**

   ```bash
   cd srcs/ubana/ubana
   git clone https://github.com/njlane314/searchingforstrangeness.git
   cd searchingforstrangeness
   source .setup
   ```

5. **Modify CMakeLists.txt:**

   ```bash
   # Add the following line to CMakeLists.txt to include searchingforstrangeness in the build 
   cd ../../ # back to src/ubana/ubana
   echo "add_subdirectory(searchingforstrangeness)" >> CMakeLists.txt
   ```

   - Integrate the cloned repository into the build process.


6. **Set up the build environment:**

   ```bash
   cd $MRB_TOP
   mrbsetenv
   ```

   - Moves to the top-level directory and sets up the build environment.

7. **Build the project:**

   ```bash
   mrb i -j4
   ```

   - Builds the project using `mrb` with 4 parallel jobs (`-j4`).

8. **Automated building**

   ```bash
   cd srcs/ubana/ubana/searchingforstrangeness
   source .container
   source .setup
   source .build
   ```
   - There exist automated bash scripts in the cloned respository that handle this configuration and build.

## Processing Files

## Configuration System

The framework uses FHiCL (Fermilab Hierarchical Configuration Language) for configuration.

### **FHiCL Structure**

### **Local Processing**

To process ROOT files locally, you can manually run a series of commands to fetch, process, combine, and clean up files. Below are the key commands and an example of how to use them.

#### Main Commands

1. **Fetch a List of Files from a SAM Dataset**:
   ```bash
   files=$(samweb list-files defname:<sam_definition> | head -n <num_files>)
   ```
   - Retrieves a specified number of files from a SAM dataset.
   Example**: `files=$(samweb list-files defname:my_dataset | head -n 5)` fetches the first 5 files.

2. **Locate the File Path**:
   ```bash
   filedir=$(samweb locate-file <file> | grep -o '/pnfs/.*' | head -n 1)
   ```
   - Finds the physical directory path of a file.
   Example**: `filedir=$(samweb locate-file somefile.root | grep -o '/pnfs/.*' | head -n 1)`.

3. **Process a File with `lar`**:
   ```bash
   lar -c <fhicl_file> -s <filepath> -T <outputfile>
   ```
   - Runs the `lar` analysis tool on a file using a FHiCL configuration.
   Example**: `lar -c analysis.fcl -s /pnfs/uboone/some/path/somefile.root -T output_1.root`.

4. **Combine Output Files**:
   ```bash
   hadd -f <combined_output> <outputfiles>
   ```
   - Merges multiple ROOT files into one.
   Example**: `hadd -f combined_output.root output_*.root`.

5. **Clean Up Temporary Files**:
   ```bash
   rm <outputfiles>
   ```
   - Deletes temporary output files.
   Example**: `rm output_*.root`.

#### Example: Process Several Files Locally
```bash
# Fetch 3 files from the dataset
files=$(samweb list-files defname:my_dataset | head -n 3)

# Process each file
counter=1
for file in $files; do
    filedir=$(samweb locate-file $file | grep -o '/pnfs/.*' | head -n 1)
    filepath="${filedir}/${file}"
    # Aftering configuring the ficl properly 
    lar -c run_eventselectionfilter.fcl -s $filepath -T output_$counter.root
    counter=$((counter + 1))
done

# Combine outputs
hadd -f combined_output.root output_*.root

# Clean up
rm output_*.root
```

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

#### Main Commands

1. **Create a Tarball**:
   ```bash
   make_tar_uboone.sh <tarball_name>
   ```
   - Packages your code into a tarball.
   Example**: `make_tar_uboone.sh my_tarball.tar`.

2. **Copy the Tarball**:
   ```bash
   cp -f <tarball_name> <tarball_dest>
   ```
   - Copies the tarball to a resilient directory.
   Example**: `cp -f my_tarball.tar /pnfs/uboone/resilient/users/myuser/`.

3. **Authenticate for Grid Access**:
   ```bash
   htgettoken -a <vault_server> -i <experiment>
   ```
   - Obtains an authentication token.
   Example**: `htgettoken -a htvaultprod.fnal.gov -i uboone`.

4. **Submit Jobs to the Grid**:
   ```bash
   project.py --xml <xml_config_file> --stage <stage> --submit
   ```
   - Submits jobs using a configuration file.
   Example**: `project.py --xml config.xml --stage analyse --submit`.

5. **Clean and Retry**:
   ```bash
   project.py --xml <xml_config_file> --stage <stage> --clean
   project.py --xml <xml_config_file> --stage <stage> --submit
   ```
   - Cleans up and retries submission.
   Example**: `project.py --xml config.xml --stage analyse --clean` followed by `project.py --xml config.xml --stage analyse --submit`.

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

`apptainer shell [options] /path/to/container`**  
    Launches an interactive shell within an Apptainer container, providing an isolated environment for running applications.  
    Example**: 
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

`source setup.sh`**  
    Executes a script to configure the shell environment, often setting variables or loading dependencies.  
    Example**: `source /setup.sh` loads necessary configurations.
    Script**: 
        ```bash
        source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh
        setup uboonecode vv10_04_07_04 -q "e26:prof"
        ```

`setup <product> <version> -q <qualifiers>`**  
    Configures the environment for a specific software product and version using UPS (Unix Product Support).  
    Example**: see above

---

### **Authentication and Grid Access**

`kx509`**  
  Converts a Kerberos ticket into an X.509 certificate for secure authentication to access grid resources.  
  Example**: Run `kx509` after `kinit` to obtain a certificate.

`voms-proxy-init -noregen -voms <virtual_organization>`**  
  Generates a proxy certificate for grid authentication tied to a virtual organization.  
  Example**: `voms-proxy-init -noregen -voms fermilab:/fermilab/uboone/Role=Analysis` creates a proxy for Fermilab.

---

### **File and Data Management**

`find /path/to/directory -type f -name "*.root"`**  
  Searches for files. 

`hadd -f output.root file1.root file2.root ...`**  
  Merges multiple ROOT files into one.  
  Example**: `hadd -f combined.root run1.root run2.root` merges two files.

`tar -cvf tarball_name.tar /path/to/code`**  
  Creates a tarball of a directory or files.  

`samweb list-files defname:<sam_definition>`**  
  Lists files in a SAM dataset.  

`samweb locate-file <file>`**  
  Finds the physical location of a SAM-managed file.  

`samweb count-files defname:<definition_name>`**  
  Counts files in a SAM dataset.  

`samweb prestage-dataset --defname=<definition_name>`**  
  Stages a SAM dataset for fast access.  

`samweb create-definition my_subset "defname:make_lambda_overlay_nohadrons_reco2_reco2 with limit 100"`**
   Creates limited dataset 

---

### **Job Submission and Monitoring**

`jobsub_q $USER`**  
  Queries the status of a user's grid jobs.  

`project.py --xml <xml_file> --stage <stage> --submit`**  
  Submits jobs to the grid, automates workflow with `*.xml` configurations.  

---

### **Data Processing and Analysis**

`lar -c <fhicl_file> -s <input_file> -T <output_file>`**  
  Runs LArSoft with a FHiCL configuration.  

`lar -c eventdump.fcl -s <your.root> | grep recob::Wire`**  
  Dumps event information from a `*.root` file and searches for `recob::Wire` objects.

---

## NuMI Input Samdefs

```bash
New_NuMI_Flux_Run_1_FHC_Pandora_Reco2_reco2_reco2      # Beam background**
prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2    # EXT/Beam-off background**
prod_strange_resample_fhc_run2_fhc_reco2_reco2         # Enriched strangeness**
make_lambda_overlay_nohadrons_reco2_reco2              # Lambda no-hadrons**
prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_v2_sample0              # Background**
cthorpe_make_k0s_events_numi_rhc_reco2_REAL_reco2_reco2                  # Kaons**
cthorpe_prod_extnumi_mcc9_v08_00_00_45_run3_run3b_reco2_all_reco2_pt1    # EXT**
cthorpe_make_hyperon_events_numi_rhc_run3b_hyperon_reco2_reco2           # Hyperons**
```

### **Splitting Input Definitions**
nl_numi_fhc_beam_run1_reco2_1000
nl_numi_fhc_ext_run1_reco2_6000
   ```bash
   samweb create-definition nl_numi_fhc_beam_run1_reco2_training_250 "defname: New_NuMI_Flux_Run_1_FHC_Pandora_Reco2_reco2_reco2 with limit 250"   
   samweb create-definition nl_numi_fhc_beam_run1_reco2_validation_250 "defname: New_NuMI_Flux_Run_1_FHC_Pandora_Reco2_reco2_reco2 with offset 250 with limit 250" 
   samweb list-files --summary "defname: nl_numi_fhc_beam_run1_reco2_validation_250"
   # The output
   File count:	250
   Total size:	427206702740
   Event count:	71153
   samweb list-files --summary "defname: nl_numi_fhc_beam_run1_reco2_training_250"
   # The output
   File count:	250
   Total size:	423907381578
   Event count:	70649
   # Even check user's list
   samweb list-definitions --user=${USER}
   ```

   ```bash
   Apptainer> samweb list-files --summary "defname: nl_ext_numi_fhc_beamoff_run1_reco2_validation_2500"
   File count:	2500
   Total size:	414497642823
   Event count:	78718
   Apptainer> samweb create-definition nl_ext_numi_fhc_beamoff_run1_reco2_training_2500 "defname: prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2 with limit 2500" 
   Dataset definition 'nl_ext_numi_fhc_beamoff_run1_reco2_training_2500' has been created with id 110354423
   Apptainer> samweb list-files --summary "defname: nl_ext_numi_fhc_beamoff_run1_reco2_training_2500"
   File count:	2500
   Total size:	367920309642
   Event count:	70654
   ```

   ```bash
   Apptainer> samweb create-definition nl_strange_numi_fhc_run2_reco2_training_982 "defname: prod_strange_resample_fhc_run2_fhc_reco2_reco2 with limit 982" 
   Dataset definition 'nl_strange_numi_fhc_run2_reco2_training_982' has been created with id 110354485
   Apptainer> samweb create-definition nl_strange_numi_fhc_run2_reco2_validation_982 "defname: prod_strange_resample_fhc_run2_fhc_reco2_reco2 with offset 982 with limit 982"
   Dataset definition 'nl_strange_numi_fhc_run2_reco2_validation_982' has been created with id 110354661
   Apptainer> samweb list-files --summary "defname: nl_strange_numi_fhc_run2_reco2_validation_982"
   File count:	982
   Total size:	564656729004
   Event count:	62935
   Apptainer> samweb list-files --summary "defname: nl_strange_numi_fhc_run2_reco2_training_982"
   File count:	982
   Total size:	1135851562244
   Event count:	126422
   ```

   ```bash
   Apptainer> samweb create-definition nl_lambda_nohadrons_reco2_validation_2000 "defname: make_lambda_overlay_nohadrons_reco2_reco2 with offset 2000 with limit 2000"
   Dataset definition 'nl_lambda_nohadrons_reco2_validation_2000' has been created with id 110354575
   Apptainer> samweb create-definition nl_lambda_nohadrons_reco2_training_2000 "defname: make_lambda_overlay_nohadrons_reco2_reco2 with limit 2000"
   Dataset definition 'nl_lambda_nohadrons_reco2_training_2000' has been created with id 110354663
   Apptainer> samweb list-files --summary "defname: nl_lambda_nohadrons_reco2_training_2000"
   File count:	2000
   Total size:	486050774712
   Event count:	71325
   Apptainer> samweb list-files --summary "defname: nl_lambda_nohadrons_reco2_validation_2000"
   File count:	2000
   Total size:	629748141141
   Event count:	91988
   ```
   
   ```bash
   [SL7][nlane@uboonegpvm03 searchingforstrangeness]$ samweb create-definition nl_prodgeni_numi_uboone_overlay_rhc_training_1600 "defname: 
   prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_v2_sample0 with limit 1600"
   Dataset definition 'nl_prodgeni_numi_uboone_overlay_rhc_training_1600' has been created with id 110403641
   [SL7][nlane@uboonegpvm03 searchingforstrangeness]$ samweb create-definition nl_prodgeni_numi_uboone_overlay_rhc_validation_1600 "defname
   : prodgenie_numi_uboone_overlay_rhc_mcc9_run3b_v28_v2_sample0 with offset 1600 with limit 1600"
   Dataset definition 'nl_prodgeni_numi_uboone_overlay_rhc_validation_1600' has been created with id 110403659
   [SL7][nlane@uboonegpvm03 searchingforstrangeness]$ samweb list-files --summary "defname:nl_prodgeni_numi_uboone_overlay_rhc_validation_1600"
   File count:     1600
   Total size:     501561617815
   Event count:    73411
   [SL7][nlane@uboonegpvm03 searchingforstrangeness]$ samweb list-files --summary "defname:nl_prodgeni_numi_uboone_overlay_rhc_training_160
   0"
   File count:     1600
   Total size:     496891362474
   Event count:    72822
   [SL7][nlane@uboonegpvm03 searchingforstrangeness]$ 
   ```



samweb list-definitions | grep -i "numi" | grep -i "run1" | grep -i "reco2" | grep -i "fhc"


EXT: nl_prod_mcc9_v08_00_00_45_extnumi_reco2_run1_all_reco2_3000