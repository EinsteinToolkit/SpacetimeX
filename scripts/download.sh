#!/bin/bash

set -ex

export SPACETIMEXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
mkdir -p "$WORKSPACE"
cd "$WORKSPACE"

# Check out Cactus
wget https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
chmod a+x GetComponents
./GetComponents --no-parallel --shallow "$SPACETIMEXSPACE/scripts/spacetimex.th"

cd Cactus

# Create a link to the SpacetimeX repository
ln -s "$SPACETIMEXSPACE" repos
# Create links for the SpacetimeX thorns
mkdir -p arrangements/SpacetimeX
pushd arrangements/SpacetimeX
ln -s ../../repos/SpacetimeX/* .
popd
