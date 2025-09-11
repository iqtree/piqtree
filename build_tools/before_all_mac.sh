# Check if running in GitHub Actions
if [ "$GITHUB_ACTIONS" = "true" ]; then
    brew update
fi

brew install llvm eigen boost libomp make

export CC=$(brew --prefix llvm)/bin/clang
export CXX=$(brew --prefix llvm)/bin/clang++
export LDFLAGS="-L$(brew --prefix libomp)/lib"
export CPPFLAGS="-I$(brew --prefix libomp)/include"
export CXXFLAGS="-I$(brew --prefix libomp)/include"

bash build_tools/build_iqtree.sh
