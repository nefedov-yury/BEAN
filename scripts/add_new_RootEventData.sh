#!/bin/bash
#
# Using:
#    add_new_RootEventData.sh  path_to_boss  RootEventData_version 

# check parameters
boss_RootEventData=`find $1/Event/RootEventData/ -maxdepth 1 -name "RootEventData*" | tail -1`
if [[ ! -d $boss_RootEventData ]]; then
  echo "check path to boss: $1"
  echo "boss_RootEventData= $boss_RootEventData"
  exit
else
  echo "copy from $boss_RootEventData"
fi

bean_RootEventData=./RootEventData_$2
if [[ -d $bean_RootEventData ]]; then
  echo "check version: $2"
  echo "dir $bean_RootEventData already exist"
  exit
else
  echo "     to $bean_RootEventData"
fi

cp -r $boss_RootEventData $bean_RootEventData

# remove unused directories and rootcint (dictionary) files:
cd $bean_RootEventData
rm -fr i386_linux26 x86_64-slc6-gcc46-opt cmt
find ./ -name CVS | xargs rm -fr
find ./ -name "RootEventData_rootcint*" | xargs rm

# copy PROOF-INF:
#cp -r ../RootEventData_6.5.5/PROOF-INF .

# create CMakeLists.txt
cat > CMakeLists.txt << EOF
INCLUDE( RootEventData )
EOF

echo " ... done"

