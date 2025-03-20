# Install dependencies using choco
choco install -y llvm eigen boost libomp make

# Set environment variables for LLVM
export CMAKE_C_COMPILER="clang"
export CMAKE_CXX_COMPILER="clang++"
export CMAKE_MAKE_PROGRAM="make"

# Build IQ-TREE
bash build_tools/build_iqtree.sh