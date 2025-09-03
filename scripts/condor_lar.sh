#! /bin/bash













































































































































































































FCL=""
INFILE=""
INLIST=""
INMODE=""
OUTFILE=""
TFILE=""
NEVT=0
NSKIP=0
SUBRUN=1
NFILE=0
NFILE_SKIP=0
NJOBS=1
NTHREADS=1
NSCHEDULES=1
ARGS=""
UPS_PRDS=""
REL=""
QUAL=""
LOCALDIR=""
LOCALTAR=""
INTERACTIVE=0
GRP=""
OUTDIR=""
LOGDIR=""
DIRSIZE=0
DIRLEVELS=0
SCRATCH=""
CLUS=""
PROC=""
PROCMAP=""
INITSCRIPT=""
INITSOURCE=""
ENDSCRIPT=""
MIDSOURCE=""
MIDSCRIPT=""
SAM_USER=$GRID_USER
SAM_GROUP=""
SAM_STATION=""
SAM_DEFNAME=""
SAM_PROJECT=""
SAM_START=0
RECUR=0
SAM_SCHEMA=""
OS=""
USE_SAM=0
MIX_DEFNAME=""
MIX_PROJECT=""
MIX_SAM=0
IFDH_OPT=""
DECLARE_IN_JOB=0
VALIDATE_IN_JOB=0
COPY_TO_FTS=0
MAINTAIN_PARENTAGE=0
EXE="lar"
INIT=""
declare -a DATAFILETYPES

while [ $
  case "$1" in


    -h|--help )
      awk '/^
      exit
      ;;


    -c|--config )
      if [ $
        FCL=$2
        shift
      fi
      ;;


    -s|--source )
      if [ $
        INFILE=$2
        shift
      fi
      ;;


    -S|--source-list )
      if [ $
        INLIST=$2
        shift
      fi
      ;;


    --inputmode )
      if [ $
        INMODE=$2
        shift
      fi
      ;;


    -o|--output )
      if [ $
        OUTFILE=$2
        shift
      fi
      ;;


    -T|--TFileName )
      if [ $
        TFILE=$2
        shift
      fi
      ;;


    -n|--nevts )
      if [ $
        NEVT=$2
        shift
      fi
      ;;


    --nskip )
      if [ $
        NSKIP=$2
        shift
      fi
      ;;


    --nfile )
      if [ $
        NFILE=$2
        shift
      fi
      ;;


    --nfile_skip )
      if [ $
        NFILE_SKIP=$2
        shift
      fi
      ;;


    --njobs )
      if [ $
        NJOBS=$2
        shift
      fi
      ;;


    --nthreads )
      if [ $
        NTHREADS=$2
        shift
      fi
      ;;


    --nschedules )
      if [ $
        NSCHEDULES=$2
        shift
      fi
      ;;


    --data_file_type )
      if [ $
        ntype=${
        DATAFILETYPES[$ntype]=$2
	shift
      fi
      ;;


    --sam_user )
      if [ $
        SAM_USER=$2
        shift
      fi
      ;;


    --sam_group )
      if [ $
        SAM_GROUP=$2
        shift
      fi
      ;;


    --sam_station )
      if [ $
        SAM_STATION=$2
        shift
      fi
      ;;


    --sam_defname )
      if [ $
        SAM_DEFNAME=$2
        USE_SAM=1
        shift
      fi
      ;;


    --sam_project )
      if [ $
        SAM_PROJECT=$2
        USE_SAM=1
        shift
      fi
      ;;


    --sam_start )
      SAM_START=1
      ;;


    --recur )
      RECUR=1
      ;;


    --sam_schema )
      if [ $
        SAM_SCHEMA=$2
        shift
      fi
      ;;


    --os )
      if [ $
        OS=$2
        shift
      fi
      ;;


    --args )
      if [ $
        shift
        ARGS=$@
        break
      fi
      ;;


    --ups )
      if [ $
        UPS_PRDS=$2
        shift
      fi
      ;;


    -r|--release )
      if [ $
        REL=$2
        shift
      fi
      ;;


    -q|-b|--build )
      if [ $
        QUAL=$2
        shift
      fi
      ;;


    --localdir )
      if [ $
        LOCALDIR=$2
        shift
      fi
      ;;


    --localtar )
      if [ $
        LOCALTAR=$2
        shift
      fi
      ;;


    --mrb )
      ;;


    --srt )
      echo "SRT run time environment is no longer supported."
      exit 1
      ;;


    -i|--interactive )
      INTERACTIVE=1
      ;;


    -g|--grid )
      ;;


    --group )
      if [ $
        GRP=$2
        shift
      fi
      ;;


    --workdir )
      if [ $
        shift
      fi
      ;;


    --outdir )
      if [ $
        OUTDIR=$2
        shift
      fi
      ;;


    --logdir )
      if [ $
        LOGDIR=$2
        shift
      fi
      ;;


    --dirsize )
      if [ $
        DIRSIZE=$2
        shift
      fi
      ;;


    --dirlevels )
      if [ $
        DIRLEVELS=$2
        shift
      fi
      ;;


    --scratch )
      if [ $
        SCRATCH=$2
        shift
      fi
      ;;


    --cluster )
      if [ $
        CLUS=$2
        shift
      fi
      ;;


    --process )
      if [ $
        PROC=$2
        shift
      fi
      ;;


    --procmap )
      if [ $
        PROCMAP=$2
        shift
      fi
      ;;


    --init-script )
      if [ $
        INITSCRIPT=$2
        shift
      fi
      ;;


    --init-source )
      if [ $
        INITSOURCE=$2
        shift
      fi
      ;;


    --end-script )
      if [ $
        ENDSCRIPT=$2
        shift
      fi
      ;;


    --mid-source )
      if [ $
        MIDSOURCE=$2
        shift
      fi
      ;;


    --mid-script )
      if [ $
        MIDSCRIPT=$2
        shift
      fi
      ;;


    --declare )
      DECLARE_IN_JOB=1
      ;;


    --validate )
      VALIDATE_IN_JOB=1
      ;;


    --copy )
      COPY_TO_FTS=1
      ;;


    --mix_defname )
      if [ $
        MIX_DEFNAME=$2
        MIX_SAM=1
        shift
      fi
      ;;


    --mix_project )
      if [ $
        MIX_PROJECT=$2
        MIX_SAM=1
        shift
      fi
      ;;


    --maintain_parentage )
      MAINTAIN_PARENTAGE=1
      ;;


    --exe )
      if [ $
        EXE=$2
        shift
      fi
      ;;


    --init )
      if [ $
        INIT=$2
        shift
      fi
      ;;


    * )
      echo "Unknown option $1"
      exit 1
  esac
  shift
done
































if [ ${
  DATAFILETYPES[0]=root
fi



echo "Nodename: `hostname -f`"
id
echo "Load average:"
cat /proc/loadavg



if [ x$QUAL = x ]; then
  QUAL="prof:e9"
fi

if [ x$SAM_GROUP = x ]; then
  SAM_GROUP=$GRP
fi

if [ x$SAM_STATION = x ]; then
  SAM_STATION=$GRP
fi



if [ x$SAM_SCHEMA = xxrootd ]; then
  SAM_SCHEMA=root
fi
if [ x$SAM_SCHEMA = xxroot ]; then
  SAM_SCHEMA=root
fi










echo "uname -r: `uname -r`"
echo "UPS_OVERRIDE: $UPS_OVERRIDE"

echo "Condor dir input: $CONDOR_DIR_INPUT"



echo "Initializing ups and mrb."

if [ x$INIT != x ]; then
  if [ ! -f $INIT ]; then
    echo "Environment initialization script $INIT not found."
    exit 1
  fi
  echo "Sourcing $INIT"
  source $INIT
else
  echo "Sourcing setup_experiment.sh"
  source ${CONDOR_DIR_INPUT}/setup_experiment.sh
fi

echo PRODUCTS=$PRODUCTS
echo "ups flavor: `ups flavor`"



unset GROUP
if [ x$GRP != x ]; then
  GROUP=$GRP
else
  echo "GROUP not specified."
  exit 1
fi
export GROUP
echo "Group: $GROUP"



echo "X509_USER_PROXY = $X509_USER_PROXY"
echo "IFDH_OPT=$IFDH_OPT"



if [ x$FCL = x ]; then
  echo "No configuration option (-c|--config) was specified."
  exit 1
fi



if [ x$OUTDIR = x ]; then
  echo "Output directory not specified."
  exit 1
fi
echo "Output directory: $OUTDIR"



if [ x$LOGDIR = x ]; then
  echo "Log directory not specified."
  exit 1
fi
echo "Log directory: $LOGDIR"






if [ $INTERACTIVE -eq 0 ]; then
  SCRATCH=$_CONDOR_SCRATCH_DIR
else
  if [ x$SCRATCH = x ]; then
    SCRATCH=$OUTDIR
  fi
fi
if [ x$SCRATCH = x -o ! -d "$SCRATCH" -o ! -w "$SCRATCH" ]; then
  echo "Local scratch directory not defined or not writable."
  exit 1
fi









TMP=`mktemp -d ${SCRATCH}/working_dir.XXXXXXXXXX`
TMP=${TMP:-${SCRATCH}/working_dir.$$}

{ [[ -n "$TMP" ]] && mkdir -p "$TMP"; } || \
  { echo "ERROR: unable to create temporary directory!" 1>&2; exit 1; }
trap "[[ -n \"$TMP\" ]] && { rm -rf \"$TMP\"; }" 0
chmod 755 $TMP
cd $TMP


echo "Scratch directory: $TMP"



echo "No longer fetching files from work directory."
echo "that's now done with using jobsub -f commands"
mkdir work
find $CONDOR_DIR_INPUT -follow -type f -exec cp {} work \;
cd work
find . -name \*.tar -exec tar xf {} \;
find . -name \*.py -exec chmod +x {} \;
find . -name \*.sh -exec chmod +x {} \;
echo "Local working directoroy:"
pwd
ls
echo



hostname > hostname.txt
echo ${CLUSTER}.${PROCESS} > jobid.txt



if [ $INTERACTIVE -ne 0 ]; then
  CLUSTER=`date +%s`
  PROCESS=0
fi



if [ x$CLUS != x ]; then
  CLUSTER=$CLUS
fi
if [ x$PROC != x ]; then
  PROCESS=$PROC
fi
if [ x$PROCMAP != x ]; then
  if [ -f $PROCMAP ]; then
    PROCESS=`sed -n $(( $PROCESS + 1 ))p $PROCMAP`
  else
    echo "Process map file $PROCMAP not found."
    exit 1
  fi
fi
if [ x$CLUSTER = x ]; then
  echo "CLUSTER not specified."
  exit 1
fi
if [ x$PROCESS = x ]; then
  echo "PROCESS not specified."
  exit 1
fi
echo "Procmap: $PROCMAP"
echo "Cluster: $CLUSTER"
echo "Process: $PROCESS"



parentdir=''
ndir=$PROCESS
while [ $DIRLEVELS -gt 0 -a $DIRSIZE -gt 0 ]; do
  parentdir=$(( $ndir % $DIRSIZE ))/$parentdir
  ndir=$(( $ndir / $DIRSIZE ))
  DIRLEVELS=$(( $DIRLEVELS - 1 ))
done
OUTPUT_SUBDIR=${parentdir}${CLUSTER}_${PROCESS}
echo "Output subdirectory: $OUTPUT_SUBDIR"



if [ ! -f $FCL ]; then
  echo "Configuration file $FCL does not exist."
  exit 1
fi



if [ x$INITSCRIPT != x ]; then
  if [ -f "$INITSCRIPT" ]; then
    chmod +x $INITSCRIPT
  else
    echo "Initialization script $INITSCRIPT does not exist."
    exit 1
  fi
fi



if [ x$INITSOURCE != x -a ! -f "$INITSOURCE" ]; then
  echo "Initialization source script $INITSOURCE does not exist."
  exit 1
fi



if [ x$ENDSCRIPT != x ]; then
  if [ -f "$ENDSCRIPT" ]; then
    chmod +x $ENDSCRIPT
  else
    echo "Finalization script $ENDSCRIPT does not exist."
    exit 1
  fi
fi



if [ x$MIDSOURCE != x -a ! -f "$MIDSOURCE" ]; then
  echo "Midstage initialization source script $MIDSOURCE does not exist."
  exit 1
fi



if [ x$MIDSCRIPT != x ]; then
  if [ -f "$MIDSCRIPT" ]; then
    chmod +x $MIDSCRIPT
  else
    echo "Midstage finalization script $MIDSCRIPT does not exist."
    exit 1
  fi
fi





if [ x$LOCALDIR != x ]; then
  mkdir $TMP/local
  cd $TMP/local



  echo "Copying local test release from directory ${LOCALDIR}."



  if [ x$IFDHC_DIR = x ]; then
    echo "Setting up ifdhc before fetching local directory."
    setup ifdhc
  fi
  echo "IFDHC_DIR=$IFDHC_DIR"
  ifdh cp -r $IFDH_OPT $LOCALDIR .
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
    exit $stat
  fi
  find . -name \*.py -exec chmod +x {} \;
  find . -name \*.sh -exec chmod +x {} \;



  cd $TMP/work
  echo "Initializing localProducts from ${LOCALDIR}."
  if [ ! -f $TMP/local/setup ]; then
    echo "Local test release directory $LOCALDIR does not contain a setup script."
    exit 1
  fi
  sed "s@setenv MRB_INSTALL.*@setenv MRB_INSTALL ${TMP}/local@" $TMP/local/setup | \
  sed "s@setenv MRB_TOP.*@setenv MRB_TOP ${TMP}@" > $TMP/local/setup.local



  if grep -q bin/shell_independence $TMP/local/setup.local; then




    echo "Setting up old version of mrb."
    unsetup mrb
    setup mrb -o
  fi



  . $TMP/local/setup.local







fi
cd $TMP/work



if [ x$LOCALTAR != x ]; then
  mkdir $TMP/local
  cd $TMP/local



  echo "Fetching test release tarball ${LOCALTAR}."



  if [ x$IFDHC_DIR = x ]; then
    echo "Setting up ifdhc before fetching tarball."
    setup ifdhc
  fi
  echo "IFDHC_DIR=$IFDHC_DIR"
  ifdh cp $LOCALTAR local.tar
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh cp failed with status ${stat}."
    exit $stat
  fi



  tar -xf local.tar



  cd $TMP/work
  echo "Initializing localProducts from tarball ${LOCALTAR}."
  sed "s@setenv MRB_INSTALL.*@setenv MRB_INSTALL ${TMP}/local@" $TMP/local/setup | \
  sed "s@setenv MRB_TOP.*@setenv MRB_TOP ${TMP}@" > $TMP/local/setup.local



  if grep -q bin/shell_independence $TMP/local/setup.local; then




    echo "Setting up old version of mrb."
    unsetup mrb
    setup mrb -o
  fi



  . $TMP/local/setup.local







fi




for prd in `echo $UPS_PRDS | tr , ' '`
do
  if ! ups active | grep -q $prd; then
    echo "Setting up $prd $REL -q ${QUAL}."
    if [ x$IFDHC_DIR != x -a x$IFBEAM_DIR = x ]; then
      unsetup ifdhc
    fi
    setup $prd $REL -q $QUAL
  fi
done

ups active

cd $TMP/work



if [ x$IFDHC_DIR = x ]; then
  echo "Setting up ifdhc again, because larsoft did not set it up."
  setup ifdhc
fi
echo "IFDH_ART_DIR=$IFDH_ART_DIR"
echo "IFDHC_DIR=$IFDHC_DIR"



if [ x$INITSCRIPT != x ]; then
  echo "Running initialization script ${INITSCRIPT}."
  ./${INITSCRIPT}
  status=$?
  if [ $status -ne 0 ]; then
    exit $status
  fi
fi

if [ x$INITSOURCE != x ]; then
  echo "Sourcing initialization source script ${INITSOURCE}."
  . $INITSOURCE
  status=$?
  if [ $status -ne 0 ]; then
    exit $status
  fi
fi



env > env.txt










rm -f condor_lar_input.list
rm -f transferred_uris.list
NFILE_TOTAL=0
parent_files=()
aunt_files=()

if [ $USE_SAM -eq 0 -a x$INFILE != x ]; then






  if [ x$INLIST != x -o $NFILE -ne 0 -o $NFILE_SKIP -ne 0 ]; then
    echo "File list options specified with single input file."
    exit 1
  fi


  parent_files=("${parent_files[@]}" $INFILE)



  NFILE_TOTAL=1
  XROOTD_URI=$INFILE
  if [ x$SAM_SCHEMA = xroot ]; then
    XROOTD_URI=`file_to_url.sh $INFILE`
  fi
  if [ $XROOTD_URI != $INFILE ]; then
    echo $INFILE > transferred_uris.list
    echo $XROOTD_URI > condor_lar_input.list
    echo "Input xrootd uri: $XROOTD_URI"
  else
    LOCAL_INFILE=`basename $INFILE`
    echo "Copying $INFILE"
    ifdh cp $INFILE $LOCAL_INFILE
    stat=$?
    if [ $stat -ne 0 ]; then
      echo "ifdh cp failed with status ${stat}."
      exit $stat
    fi
    if [ -f $LOCAL_INFILE -a $stat -eq 0 ]; then
      echo $INFILE > transferred_uris.list
      echo $LOCAL_INFILE > condor_lar_input.list
    else
      echo "Error fetching input file ${INFILE}."
      exit 1
    fi
  fi

elif [ $USE_SAM -eq 0 -a x$INLIST != x ]; then





  if [ ! -f $INLIST ]; then
    echo "Input file list $INLIST does not exist."
    exit 1
  fi



  NFILE_TOTAL=`cat $INLIST | wc -l`
  echo "Input file list contains $NFILE_TOTAL total files."





  MAX_TOTAL=$(( $NFILE * $NJOBS ))
  if [ $MAX_TOTAL -gt 0 -a $NFILE_TOTAL -gt $MAX_TOTAL ]; then
    NFILE_TOTAL=$MAX_TOTAL
    echo "Number of files to be processed will be limited to ${NFILE_TOTAL}."
  fi




  if [ $NJOBS -ne 0 ]; then



    if [ $NFILE_SKIP -ne 0 ]; then
      echo "Illegal options specified with --njobs."
      exit 1
    fi





    MYNJOBS=$NJOBS
    if [ $MYNJOBS -gt $NFILE_TOTAL ]; then
      MYNJOBS=$NFILE_TOTAL
    fi



    NFILE_SKIP=$(( $PROCESS * $NFILE_TOTAL / $MYNJOBS ))
    MYNFILE=$(( ( $PROCESS + 1 ) * $NFILE_TOTAL / $MYNJOBS - $NFILE_SKIP ))
    if [ $MYNFILE -eq 0 -o $NFILE_SKIP -ge $NFILE_TOTAL ]; then
      echo "This worker did not get any input files."
      exit 1
    fi
    if [ $MYNFILE -lt $NFILE -o $NFILE -eq 0 ]; then
      NFILE=$MYNFILE
    fi
  fi



  echo "Skipping $NFILE_SKIP files."
  if [ $NFILE -eq 0 ]; then
    echo "Processing all remaining files."
  else
    echo "Processing $NFILE files."
  fi



  nfile=0
  nfskip=$NFILE_SKIP
  nmax=$NFILE
  while read infile; do
    if [ $nfskip -gt 0 ]; then
      nfskip=$(( $nfskip - 1 ))
    else




      if [ ! -f condor_lar_input.list ]; then
        touch condor_lar_input.list
      fi

      XROOTD_URI=$infile
      if [ x$SAM_SCHEMA = xroot ]; then
        XROOTD_URI=`file_to_url.sh $infile`
      fi
      if [ $XROOTD_URI != $infile ]; then
        echo $infile >> transferred_uris.list
        echo $XROOTD_URI >> condor_lar_input.list
        echo "Input xrootd uri: $XROOTD_URI"
      else
        LOCAL_INFILE=`basename $infile`
        if grep -q $LOCAL_INFILE condor_lar_input.list; then
          LOCAL_INFILE=input${nfile}.root
	  if [ "$INMODE" = 'textfile' ]; then
	    LOCAL_INFILE=input${nfile}.txt
	  fi
        fi
        echo "Copying $infile"
        ifdh cp $infile $LOCAL_INFILE
        stat=$?
        if [ $stat -ne 0 ]; then
          echo "ifdh cp failed with status ${stat}."
          exit $stat
        fi
        if [ -f $LOCAL_INFILE -a $stat -eq 0 ]; then
          echo $infile >> transferred_uris.list
          echo $LOCAL_INFILE >> condor_lar_input.list
	  parent_files=("${parent_files[@]}" $LOCAL_INFILE)
        else
          echo "Error fetching input file ${infile}."
          exit 1
        fi
      fi
      nmax=$(( $nmax - 1 ))
      if [ $nmax -eq 0 ]; then
        break
      fi
    fi
    nfile=$(( $nfile + 1 ))
  done < $INLIST
fi

NFILE_LOCAL=0
if [ $USE_SAM -eq 0 -a x$SAM_SCHEMA != xroot ]; then
  if [ -f condor_lar_input.list ]; then





    xargs ls -s1 < condor_lar_input.list | sort -nr | awk '{print $2}' > newcondor_lar_input.list
    mv -f newcondor_lar_input.list condor_lar_input.list
    echo "Local input file list:"
    cat condor_lar_input.list
    NFILE_LOCAL=`cat condor_lar_input.list | wc -l`
  else
    echo "No local input files."
  fi
  echo "Local input list has $NFILE_LOCAL files."
fi


nfcls=0

while read -r line
do

  if [ "$(echo $line | awk '{print $1}')" = "
    stage="$(echo $line | awk '{print $2}')"
    stage_fcl="Stage$stage.fcl"
    nfcls=$(( $nfcls + 1 ))
    continue
  fi

  if [ "$line" = "

    continue
  fi
  echo $line >> $stage_fcl
done < $FCL


stage=0

echo "Start loop over stages"
while [ $stage -lt $nfcls ]; do
  FCL="Stage$stage.fcl"







  if [ $stage -eq 0 -a $USE_SAM -eq 0 ] && [ $NFILE_TOTAL -eq 0 -o "$INMODE" = 'textfile' ]; then




    if [ $NSKIP -gt 0 ]; then
      echo "Illegal option --nskip specified with no input."
      exit 1
    fi



    NSKIP=$(( $PROCESS * $NEVT / $NJOBS ))
    NEV=$(( ( $PROCESS + 1 ) * $NEVT / $NJOBS - $NSKIP ))
    NSKIP=0
    NEVT=$NEV



    SUBRUN=$(( $PROCESS + 1))
    cat <<EOF > subrun_wrapper.fcl


source.firstSubRun: $SUBRUN

EOF
    if [ "$INMODE" = 'textfile' ]; then

      if [ $NFILE_LOCAL -ne 1 ]; then
        echo "Text file input mode specified with wrong number of input files."
        exit 1
      fi
      echo "physics.producers.generator.InputFileName: \"`cat condor_lar_input.list`\"" >> subrun_wrapper.fcl
    fi

    FCL=subrun_wrapper.fcl

    echo "MC subrun: $SUBRUN"
    echo "Number of MC events: $NEVT"

  fi



  PURL=''
  CPID=''
  if [ $USE_SAM -ne 0 -a $stage -eq 0 ]; then
    echo "In SAM if"



    if [ x$SAM_PROJECT = x ]; then
      echo "No sam project was specified."
      exit 1
    fi
    echo "Sam project: $SAM_PROJECT"



    if [ $SAM_START -ne 0 ]; then
      if [ x$SAM_DEFNAME != x ]; then






        nf=`ifdh translateConstraints "defname: $SAM_DEFNAME" | wc -l`
        if [ $nf -eq 0 ]; then
          echo "Input dataset $SAM_DEFNAME is empty."
          exit 1
        fi
        if [ $NFILE -ne 0 -a $nf -gt $NFILE ]; then
          limitdef=${SAM_PROJECT}_limit_$NFILE




          existdef=`ifdh describeDefinition $limitdef 2>/dev/null | grep 'Definition Name:' | wc -l`
          if [ $existdef -gt 0 ]; then
            echo "Using already created limited dataset definition ${limitdef}."
          else
            ifdh createDefinition $limitdef "defname: $SAM_DEFNAME with limit $NFILE" $SAM_USER $SAM_GROUP



            echo "Created limited dataset definition ${limitdef}."
          fi




          SAM_DEFNAME=$limitdef
        fi



          if [ $RECUR -ne 0 ]; then
            echo "Forcing snapshot"
            SAM_DEFNAME=${SAM_DEFNAME}:force
          fi



        echo "Starting project $SAM_PROJECT using sam dataset definition $SAM_DEFNAME"
        ifdh startProject $SAM_PROJECT $SAM_STATION $SAM_DEFNAME $SAM_USER $SAM_GROUP
        if [ $? -eq 0 ]; then
          echo "Start project succeeded."
        else
          echo "Start projet failed."
          exit 1
        fi
      fi

      if [ x$SAM_DEFNAME = x ]; then

        echo "Start project requested, but no definition was specified."
        exit 1
      fi

    fi






    PURL=`ifdh findProject $SAM_PROJECT $SAM_STATION`
    if [ x$PURL = x ]; then
      echo "Unable to find url for project ${SAM_PROJECT}."
      exit 1
    else
      echo "Project url: $PURL"
    fi



    NODE=`hostname`
    APPFAMILY=art




    APPNAME=`fhicl-dump $FCL | grep process_name: | head -1 | tr -d '"' | awk '{print $2}'`
    if [ $? -ne 0 ]; then
      echo "fhicl-dump $FCL failed to run. May be missing a ups product, library, or fcl file."
      exit 1
    fi
    if [ x$APPNAME = x ]; then
      echo "Trouble determining application name."
      echo "cat $FCL"
      cat $FCL
      exit 1
    fi



    if [ x$REL = x ]; then
      REL=1
    fi




    DESC=$JOBSUBJOBID
    if [ x$DESC = x ]; then
      DESC=$FCL
    fi

    echo "Starting consumer process."
    echo "ifdh establishProcess $PURL $APPNAME $REL $NODE $SAM_USER $APPFAMILY $DESC $NFILE $SAM_SCHEMA"
    CPID=`ifdh establishProcess $PURL $APPNAME $REL $NODE $SAM_USER $APPFAMILY $DESC $NFILE $SAM_SCHEMA`
    if [ x$CPID = x ]; then
      echo "Unable to start consumer process for project url ${PURL}."
      exit 1
    else
      echo "Consumer process id $CPID"
    fi




    echo $SAM_PROJECT > sam_project.txt
    echo $CPID > cpid.txt

  fi



  if [ $MIX_SAM -ne 0 ]; then
    echo "In Mix SAM if"



    if [ x$MIX_PROJECT = x ]; then
      echo "No mix sam project was specified."
      exit 1
    fi
    echo "Mix project: $MIX_PROJECT"



    if [ $SAM_START -ne 0 ]; then
      if [ x$MIX_DEFNAME != x ]; then

        echo "Starting project $MIX_PROJECT using sam dataset definition $MIX_DEFNAME"
        ifdh startProject $MIX_PROJECT $SAM_STATION $MIX_DEFNAME $SAM_USER $SAM_GROUP
        if [ $? -eq 0 ]; then
          echo "Start project succeeded."
        else
          echo "Start projet failed."
          exit 1
        fi
      fi

      if [ x$MIX_DEFNAME = x ]; then

        echo "Start project requested, but no mix definition was specified."
        exit 1
      fi
    fi
  fi







  LAROPT="-c $FCL --rethrow-default"
  echo "Laropt: $LAROPT"
  if [ -f condor_lar_input.list -a $stage -eq 0 ]; then
    if [ "$INMODE" != 'textfile' ]; then
      LAROPT="$LAROPT -S condor_lar_input.list"

    fi
  fi



  if echo $OUTFILE | grep -q :; then
    outfile=''
  else
    outfile=$OUTFILE
  fi
  field=$(( $stage + 1 ))
  outfile_stage=`echo $OUTFILE | cut -d: -f$field`
  if [ x$outfile_stage != x ]; then
    outfile=$outfile_stage
  fi
  if [ x$outfile != x ]; then
    LAROPT="$LAROPT -o `basename $outfile .root`$stage.root"
    outstem=`basename $OUTFILE .root`
  fi

  if [ x$TFILE != x ]; then
    LAROPT="$LAROPT -T $TFILE"
  fi

  if [ $NEVT -ne 0 ]; then
    LAROPT="$LAROPT -n $NEVT"
  fi

  if [ $NSKIP -ne 0 ]; then
    LAROPT="$LAROPT --nskip $NSKIP"
  fi

  if [ $NTHREADS -gt 1 ]; then
    LAROPT="$LAROPT --nthreads $NTHREADS"
  fi

  if [ $NSCHEDULES -gt 1 ]; then
    LAROPT="$LAROPT --nschedules $NSCHEDULES"
  fi

  if [ x$PURL != x -a $stage -eq 0 ]; then
    LAROPT="$LAROPT --sam-web-uri $PURL"
  fi

  if [ x$CPID != x -a $stage -eq 0 ]; then
    LAROPT="$LAROPT --sam-process-id $CPID"
  fi

  if [ -n "$ARGS" ]; then
    LAROPT="$LAROPT $ARGS"
  fi



  if [ x$MIDSOURCE != x ]; then
    echo "Sourcing midstage initialization source script ${MIDSOURCE}."
    . $MIDSOURCE
    status=$?
    if [ $status -ne 0 ]; then
      exit $status
    fi
  fi

  if [ $stage -ne 0 ]; then
    LAROPT="$LAROPT -s $next_stage_input"
  fi



  env > env${stage}.txt



  fhicl-dump $FCL > cfgStage$stage.fcl



  echo
  echo "Proxy:"
  echo
  voms-proxy-info -all


  pwd



  if echo $EXE | grep -q :; then
    exe='lar'
  else
    exe=$EXE
  fi
  field=$(( $stage + 1 ))
  exe_stage=`echo $EXE | cut -d: -f$field`
  if [ x$exe_stage != x ]; then
    exe=$exe_stage
  fi
  echo "$exe $LAROPT"
  echo "$exe $LAROPT" > commandStage$stage.txt
  $exe $LAROPT > larStage$stage.out 2> larStage$stage.err
  stat=$?
  echo $stat > larStage$stage.stat
  echo "$exe completed with exit status ${stat}."
  if [ $stat -ne 0 ]; then
    echo
    echo "Proxy:"
    echo
    voms-proxy-info -all
    echo
    echo "tail -1000 larStage$stage.out"
    echo
    tail -1000 larStage$stage.out
    echo
    echo "tail -1000 larStage$stage.err"
    echo
    tail -1000 larStage$stage.err
    echo
  fi



  if [ $USE_SAM -ne 0 -a $stage -eq 0 ]; then



    if [ x$CPID = x -a -f cpid.txt ]; then
      CPID=`cat cpid.txt`
    fi
    ifdh translateConstraints "consumer_process_id $CPID and consumed_status consumed" > consumed_files.list



    ifdh endProcess $PURL $CPID



    nprj=`ifdh translateConstraints "snapshot_for_project_name $SAM_PROJECT" | wc -l`
    nconsumed=`ifdh translateConstraints "project_name $SAM_PROJECT and consumed_status consumed" | wc -l`
    echo "$nprj files in project, $nconsumed files consumed so far."

    if [ $SAM_START -ne 0 -o \( $nprj -gt 0 -a $nconsumed -eq $nprj \) ]; then
      echo "Stopping project."
      ifdh endProject $PURL
    fi
  fi


  if [ $stat -ne 0 ]; then
    break
  fi



  if [ x$MIDSCRIPT != x ]; then
    echo "Running midstage finalization script ${MIDSCRIPT}."
    ./${MIDSCRIPT} $stage
    status=$?
    if [ $status -ne 0 ]; then
      exit $status
    fi
  fi



  if [ $stage -ne 0 ]; then
   rm -rf $next_stage_input
  fi






  next_stage_input=`ls -t1 *.root | egrep -v 'celltree|hist|larlite|larcv|Supplemental|TGraphs' | artroot_filter.py | head -n1`



  nc=`echo $next_stage_input | wc -c`
  if [ $nc -ge 200 ]; then
    base=`basename $next_stage_input`
    ext=${base
    stem=${base%.*}
    newstem=`echo $stem | cut -c1-150`_`uuidgen`
    echo "mv $next_stage_input ${newstem}.${ext}"
    mv $next_stage_input ${newstem}.${ext}
    next_stage_input=${newstem}.${ext}
  fi

  mixed_files=`sam_metadata_dumper $next_stage_input | grep mixparent | awk -F ":" '{gsub("\"" ,""); gsub(",",""); gsub(" ",""); print $2}' | sort -u`

  if [ x"$mixed_files" != x ]; then
    aunt_files=("${aunt_files[@]}" $mixed_files)
  fi

  stage=$[$stage +1]



  if [ -f time.db ]; then
    mv time.db time$stage.db
  fi
  if [ -f mem.db ]; then
    mv mem.db mem$stage.db
  fi

done





if [ $MIX_SAM -ne 0 ]; then



  if [ $SAM_START -ne 0 ]; then
    echo "Stopping project."
    MURL=`ifdh findProject $MIX_PROJECT $SAM_STATION`
    ifdh endProject $MURL
  fi
fi



if [ $USE_SAM -eq 0 -a x$SAM_SCHEMA != xroot -a -f condor_lar_input.list ]; then
  while read file; do
    rm -f $file
  done < condor_lar_input.list
fi



if [ x$ENDSCRIPT != x ]; then
  echo "Running end-of-job script ${ENDSCRIPT}."
  ./${ENDSCRIPT}
  status=$?
  if [ $status -ne 0 ]; then
    exit $status
  fi
fi















ran=0
if [ $USE_SAM -eq 0 -a x$INFILE = x -a x$INLIST = x ]; then
  ran=1
fi

for ftype in ${DATAFILETYPES[*]}; do
  for datafile in *.${ftype}; do
    if [ -f $datafile ]; then
      nc=`echo $datafile | wc -c`
      if ifdh getMetadata $datafile > /dev/null 2> /dev/null; then
        nc=200
      fi
      if [ -f ${datafile}.json -o $ran != 0 -o $nc -ge 200 ]; then
        base=`basename $datafile`
        ext=${base
        stem=${base%.*}
        newstem=`echo $stem | cut -c1-150`_`uuidgen`
        echo "mv $datafile ${newstem}.${ext}"
        mv $datafile ${newstem}.${ext}
        if [ -f ${datafile}.json ]; then
          mv ${datafile}.json ${newstem}.${ext}.json
        fi
      fi
    fi
  done
done




for ftype in ${DATAFILETYPES[*]}; do
  for datafile in *.${ftype}; do
    if [ -f $datafile ]; then
      json=${datafile}.json
      if [ -f $json ]; then
        ./root_metadata.py --output="${json}2" "$datafile" >& /dev/null
        ./merge_json.py $json ${json}2 > ${json}3
        mv -f ${json}3 $json
        rm ${json}2
      else
        ./root_metadata.py --output="$json" "$datafile" >& /dev/null
      fi
    fi
  done
done


stageStat=0
overallStat=0
while [ $stageStat -lt $nfcls ]; do
  stat=`cat larStage$stageStat.stat`
  if [[ "$stat" = 65 && $ART_VERSION < v2_01 ]]; then

    for json in *.json; do
        if grep -q '"events": *"0"' $json; then
         stat=0
        fi
    done
  fi
  overallStat=$[$stat+$overallStat]



  stageStat=$[$stageStat +1]
done
echo $overallStat > lar.stat
valstat=$overallStat



mkdir out
mkdir log





for ftype in ${DATAFILETYPES[*]}; do
  for datafile in *.${ftype}; do
    if [ -f $datafile ]; then
      mv $datafile out
      if [ -f ${datafile}.json ]; then
        mv ${datafile}.json log
      fi
    fi
  done
done



for outfile in *; do
  if [ -f $outfile ]; then
    mv $outfile log
  fi
done



if [ $VALIDATE_IN_JOB -eq 1 ]; then

    if [ $USE_SAM -ne 0 ]; then
      id=`cat log/cpid.txt`
      parent_files=($(ifdh translateConstraints "consumer_process_id=$id and consumed_status consumed"))
      stat=$?
      if [ $stat -ne 0 ]; then
        echo "Failed to determine parentage."

      fi
    fi

    echo "The file's parents are: "

    for elt in ${parent_files[*]};
    do
      echo $elt
    done

    echo "The file's aunts are: "
    for elt in ${aunt_files[*]};
    do
      echo $elt
    done



    if [ $MAINTAIN_PARENTAGE -eq 1 ]; then
       export JOBS_PARENTS=`echo ${parent_files[*]}`
       export JOBS_AUNTS=`echo ${aunt_files[*]}`
    fi



    valstat=$overallStat
    if [ $valstat -eq 0 ]; then
      curdir=`pwd`
      cd $curdir/log
      dataopt=''
      for ftype in ${DATAFILETYPES[*]}; do
        dataopt="$dataopt --data_file_type $ftype"
      done
      echo "./validate_in_job.py --dir $curdir/out --logfiledir $curdir/log --outdir $OUTDIR/$OUTPUT_SUBDIR --declare $DECLARE_IN_JOB --copy $COPY_TO_FTS --maintain_parentage $MAINTAIN_PARENTAGE $dataopt"
      ./validate_in_job.py --dir $curdir/out --logfiledir $curdir/log --outdir $OUTDIR/$OUTPUT_SUBDIR --declare $DECLARE_IN_JOB --copy $COPY_TO_FTS --maintain_parentage $MAINTAIN_PARENTAGE $dataopt
      valstat=$?
      cd $curdir
    fi

fi



rm -f log.tar
tar -cjf log.tar -C log .
mv log.tar log




echo "Setting up current version of ifdhc."
if [ x$IFDHC_DIR != x ]; then
  unsetup ifdhc
fi
setup ifdhc
echo "IFDHC_DIR=$IFDHC_DIR"



export IFDH_CP_MAXRETRIES=5

echo "Make directory ${LOGDIR}/${OUTPUT_SUBDIR}."
date
echo "ifdh mkdir_p $IFDH_OPT ${LOGDIR}/$OUTPUT_SUBDIR"
ifdh mkdir_p $IFDH_OPT ${LOGDIR}/$OUTPUT_SUBDIR 2> err.txt
stat=$?
if [ $stat -ne 0 ]; then
  echo "ifdh command returned status $stat"
  cat err.txt
fi
rm -f err.txt
echo "Done making directory ${LOGDIR}/${OUTPUT_SUBDIR}."
date

if [ ${OUTDIR} != ${LOGDIR} ]; then
  echo "Make directory ${OUTDIR}/${OUTPUT_SUBDIR}."
  date
  echo "ifdh mkdir_p $IFDH_OPT ${OUTDIR}/$OUTPUT_SUBDIR"
  ifdh mkdir_p $IFDH_OPT ${OUTDIR}/$OUTPUT_SUBDIR 2> err.txt
  stat=$?
  if [ $stat -ne 0 ]; then
    echo "ifdh command returned status $stat"
    cat err.txt
  fi
  rm -f err.txt
  echo "Done making directory ${OUTDIR}/${OUTPUT_SUBDIR}."
  date
fi



statout=0
echo "ls log"
ls log
echo "ifdh cp -D $IFDH_OPT log/log.tar ${LOGDIR}/$OUTPUT_SUBDIR"
ifdh cp -D $IFDH_OPT log/log.tar ${LOGDIR}/$OUTPUT_SUBDIR
date
stat=$?
if [ $stat -ne 0 ]; then
  statout=1
  echo "ifdh cp failed with status ${stat}."
fi



if [ $COPY_TO_FTS -eq 0 ]; then

  if [ "$( ls -A out )" ]; then
    echo "ifdh cp -D $IFDH_OPT out/* ${OUTDIR}/$OUTPUT_SUBDIR"
    ifdh cp -D $IFDH_OPT out/* ${OUTDIR}/$OUTPUT_SUBDIR
    stat=$?
    if [ $stat -ne 0 ]; then
      statout=1
      echo "ifdh cp failed with status ${stat}."
    fi
  fi

fi

if [ $statout -eq 0 -a -f log/lar.stat ]; then
  statout=`cat log/lar.stat`
fi

if [ $statout -eq 0 ]; then
  statout=$valstat
fi

exit $statout