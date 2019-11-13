#!/usr/bin/env bash

# Supported genomes
GENOMES=("hsapi19" "hsapi38" "mmusc10")

# Help function
function print_usage() {
  echo "USAGE: $0 tool databasename data" >&2
  echo "  tool: build, run, manual" >&2
  echo "    build: Build RNA database into a docker volume" >&2
  echo "    run: Run TRAX with a prebuilt docker volume, must also supply the following arguments" >&2
  echo "      experiment_name: The name for your samples that will be made as a folder in the data parameter" >&2
  echo "      sample_file: The supplied sample file located under the data parameter" >&2
  echo "    manual: Run container with prebuilt docker volume" >&2
  echo "  databasename: hsapi19, hsapi38, mmusc10" >&2
  echo "  data: Directory to mount with the data (optional)" >&2 
}

# Function to start trax automatically
function docker_trax() {
  if [ -z "$4" ]
  then
    print_usage
  else
    docker run --rm -it --name trax-${USER} \
      --user=`id -u`:`id -g` \
      -v rnadb-${1}:/rnadb \
      -v ${2}:/data \
      -v `pwd`/output:/trax/output \
      -v `pwd`/sample:/trax/samples \
      trax ./processsamples.py \
      --experimentname=/output \
      --databasename=/rnadb/db \
      --samplefile=/data/${4} \
      --ensemblgtf=/rnadb/genes.gtf
  fi
}

# Function to start trax container
function docker_manual() {
  docker run --rm -it --name trax-${USER} \
    --user=`id -u`:`id -g` \
    -v rnadb-${1}:/rnadb \
    -v ${2}:/data \
    trax
}

# Function to start the container with a named volume
function docker_build_db() {
  docker volume create rnadb-${1}
  docker run --rm -it --name trax-build-rnadb-${1} \
    -v rnadb-${1}:/rnadb \
    -v ${2}:/data \
    trax ./quickdb.bash ${1}
}

# Init function
if test ${1} = "build"
then
  [[ ${GENOMES[*]} =~ ${2} ]] && docker_build_db ${2} ${3} || print_usage
elif test ${1} = "run"
then
  [[ ${GENOMES[*]} =~ ${2} ]] && docker_trax ${2} ${3} ${4} ${5} || print_usage
elif test ${1} = "manual"
then
  [[ ${GENOMES[*]} =~ ${2} ]] && docker_manual ${2} ${3} || print_usage
else
  print_usage
fi

