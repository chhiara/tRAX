#!/usr/bin/env bash

# Supported genomes
GENOMES=("hsapi19" "hsapi38" "mmusc10")

# Help function
function print_usage() {
  echo "USAGE: $0 tool databasename data" >&2
  echo "  tool: make, build, run, manual" >&2
  echo "    make: Build the docker container" >&2
  echo "    build: Build RNA database into a docker volume" >&2
  echo "    run: Run TRAX with a prebuilt docker volume, must also supply the following arguments" >&2
  echo "      experiment_name: The name for your samples that will be made as a folder in the data parameter" >&2
  echo "      sample_file: The supplied sample file located under the data parameter" >&2
  echo "    manual: Run container with prebuilt docker volume" >&2
  echo "  databasename: hsapi19, hsapi38, mmusc10" >&2
  echo "  data: Directory to mount with the data (optional)" >&2 
}

# Function to start trax automatically and run
#function docker_trax() {
#  if [ -z "$4" ]
#  then
#    print_usage
#  else
#    docker run --rm -it --name trax-${USER} \
#      --user=`id -u`:`id -g` \
#      -v rnadb-${1}:/rnadb \
#      -v ${2}:${2} \
#      trax ./processsamples.py \
#      --experimentname=/home/jerry/data/output \
#      --databasename=/rnadb/db \
#      --samplefile=/home/jerry/data/${4} \
#      --ensemblgtf=/rnadb/genes.gtf
#  fi
#}

# Function to build the Docker container
function docker_make() {
  docker build -f Dockerfile -t trax .
}

# Function to start the container and build a RNA database
function docker_build_db() {
  docker volume create rnadb-${1}
  docker run --rm -it --name trax-build-rnadb-${1} \
    -v rnadb-${1}:/rnadb \
    trax ./quickdb_trax.bash ${1}
}

# Function to start a manual Docker TRAX container
function docker_manual() {
  docker run --rm -it --name trax-${USER}-2 \
  --user=`id -u`:`id -g` \
  -v rnadb-${1}:/rnadb \
  #for i in ${@:2}
  #do
  #  -v  $i:$i
  #done
  trax
}

# Init function
if test ${1} = "make"
then
 docker_make
elif test ${1} = "build"
then
  [[ ${GENOMES[*]} =~ ${2} ]] && docker_build_db ${@:2} || print_usage
#elif test ${1} = "run"
#then
#  [[ ${GENOMES[*]} =~ ${2} ]] && docker_trax ${@:2} || print_usage
elif test ${1} = "manual"
then
  [[ ${GENOMES[*]} =~ ${2} ]] && docker_manual ${@:2} || print_usage
else
  print_usage
fi

