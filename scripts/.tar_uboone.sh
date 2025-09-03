#! /bin/bash




































function dohelp {
  echo "Usage: make_tar_uboone.sh [-h|--help] [-d dir] <tarball-name>"
}



if [ $
  dohelp
  exit
fi



tar=''
dir=''
if [ x$SRT_PRIVATE_CONTEXT != x ]; then
  dir=$SRT_PRIVATE_CONTEXT
fi
if [ x$MRB_INSTALL != x ]; then
  dir=$MRB_INSTALL
fi
if [ x$dir = x ]; then
  dir=`pwd`
fi

while [ $
  case "$1" in


    -h|--help )
      dohelp
      exit
      ;;


    -d )
      if [ $
        dir=$2
        shift
      fi
      ;;


    -* )
      echo "Unrecognized option $1"
      dohelp
      exit
      ;;


    * )
      if [ x$tar = x ]; then
        tar=$1
      else
        echo "Too many arguments."
        dohelp
        exit 1
      fi

  esac
  shift
done



if [ x$dir = x ]; then
  echo "No source directory specified."
  dohelp
  exit 1
fi
if [ ! -d $dir ]; then
  echo "Directory $dir doesn't exist."
  exit 1
fi



if [ x$tar = x ]; then
  echo "No tarball specified."
  dohelp
  exit 1
fi



if [ -f $tar ]; then
  rm $tar
fi

ls -A $dir | egrep -v '.root$|tmp' | tar -C $dir -T- -czf $tar --exclude=.svn --exclude=.git --exclude=env --exclude=scratch --exclude=\*.tar



