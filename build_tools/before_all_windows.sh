# Install dependencies using choco

export Boost_INCLUDE_DIR=$(echo $Boost_INCLUDE_DIR | sed 's|\\|/|g')
export Boost_LIBRARY_DIRS=$(echo $Boost_LIBRARY_DIRS | sed 's|\\|/|g')

echo "Boost_INCLUDE_DIR: $Boost_INCLUDE_DIR"
echo "Boost_LIBRARY_DIRS: $Boost_LIBRARY_DIRS"

choco install -y llvm --version=14.0.6 --allow-downgrade
choco install -y eigen 

llvm_tools_dir=("${ChocolateyInstall}/lib/llvm/tools/"*/)
LLVM_BIN="${llvm_tools_dir%/}/bin"

export PATH="$LLVM_BIN:$PATH"

echo "Using LLVM from: $LLVM_BIN"

# Build IQ-TREE
bash build_tools/build_iqtree.sh