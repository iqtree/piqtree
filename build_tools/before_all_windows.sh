echo "BOOST_ROOT=$BOOST_ROOT"
echo "Boost_INCLUDE_DIR=$Boost_INCLUDE_DIR"
echo "Boost_LIBRARY_DIRS=$Boost_LIBRARY_DIRS"
echo "CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH"

# Install dependencies using choco
choco install -y llvm eigen make

# Set environment variables for LLVM
export CMAKE_C_COMPILER="clang"
export CMAKE_CXX_COMPILER="clang++"
export CMAKE_MAKE_PROGRAM="make"

# Build IQ-TREE
bash build_tools/build_iqtree.sh