# Install dependencies using choco

export Boost_INCLUDE_DIR=$(echo $Boost_INCLUDE_DIR | sed 's|\\|/|g')
export Boost_LIBRARY_DIRS=$(echo $Boost_LIBRARY_DIRS | sed 's|\\|/|g')

echo "Boost_INCLUDE_DIR: $Boost_INCLUDE_DIR"
echo "Boost_LIBRARY_DIRS: $Boost_LIBRARY_DIRS"

choco install -y llvm --version=14.0.6 --allow-downgrade
choco install -y eigen 

CHOCO_UNIX="$(cygpath "${ChocolateyInstall}")"
echo "Chocolatey install root: $CHOCO_UNIX"

echo "Listing contents of ${CHOCO_UNIX}/lib/llvm/tools:"
find "${CHOCO_UNIX}/lib/llvm/tools" || true

# Build IQ-TREE
bash build_tools/build_iqtree.sh
