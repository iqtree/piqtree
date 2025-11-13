# Check if running in GitHub Actions
if [ "$GITHUB_ACTIONS" = "true" ]; then
    brew update
fi

brew install llvm eigen@3 boost libomp make

export LDFLAGS="-L$(brew --prefix libomp)/lib"
export CPPFLAGS="-I$(brew --prefix libomp)/include -I$(brew --prefix eigen@3)/include"
export CXXFLAGS="-I$(brew --prefix libomp)/include"

bash build_tools/build_iqtree.sh
