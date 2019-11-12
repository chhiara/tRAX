#!/usr/bin/env bash

# Supported genomes
GENOMES=("hsapi19" "hsapi38" "mmusc10")

# Help function
function print_usage() {
  echo "USAGE: $0 tool databasename" >&2
  echo "  tool: builddb, runcontainer" >&2
  echo "  databasename: hsapi19, hsapi38, mmusc10" >&2
}

function docker_trax() {
  docker run --rm -it --name trax-${USER} \
    -v rnadb-${2}:/rnadb \
    -v `pwd`:/trax \
    trax
}

# Function to start the container with a named volume
function docker_builddb() {
  docker volume create rnadb-${2}
  docker run --rm -it --name trax-build-rnadb-${2} \
    -v rnadb-${2}:/rnadb \
    -v `pwd`:/trax \
    trax ./quickdb.bash ${2}
}

# Init function
if test "${1}" = "builddb"
then
  [[ $GENOMES =~ (^|[[:space:]])${2}($|[[:space:]]) ]] && docker_builddb || print_usage
elif test "${1}" = "runcontainer"
then
  [[ $GENOMES =~ (^|[[:space:]])${2}($|[[:space:]]) ]] && docker_trax || print_usage
else
  print_usage
fi

