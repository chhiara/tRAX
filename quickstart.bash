#!/usr/bin/env bash

# Supported genomes
GENOMES=("hsapi19" "hsapi38" "mmusc10")

# Help function
function print_usage() {
  echo "USAGE: $0 tool databasename" >&2
  echo "  tool: build, run, manual" >&2
  echo "    build: Build RNA database into a docker volume" >&2
  echo "    run: Run TRAX with a prebuilt docker volume" >&2
  echo "    manual: Run container with prebuilt docker volume" >&2
  echo "  databasename: hsapi19, hsapi38, mmusc10" >&2
}

# Function to start trax automatically
function docker_trax() {
  docker run --rm -it --name trax-${USER} \
    -v rnadb-${1}:/rnadb \
    -v `pwd`:/trax \
    trax ./processsamples.py
}

# Function to start trax container
function docker_manual() {
  docker run --rm -it --name trax-${USER} \
    -v rnadb-${1}:/rnadb \
    -v `pwd`:/trax \
    trax
}

# Function to start the container with a named volume
function docker_build_db() {
  docker volume create rnadb-${1}
  docker run --rm -it --name trax-build-rnadb-${1} \
    -v rnadb-${1}:/rnadb \
    -v `pwd`:/trax \
    trax ./quickdb.bash ${1}
}

# Init function
if test ${1} = "build"
then
  [[ ${GENOMES[*]} =~ ${2} ]] && docker_build_db ${2} || print_usage
elif test ${1} = "run"
then
  [[ ${GENOMES[*]} =~ ${2} ]] && docker_trax ${2} || print_usage
elif test ${1} = "manual"
then
  [[ ${GENOMES[*]} =~ ${2} ]] && docker_manual ${2} || print_usage
else
  print_usage
fi

