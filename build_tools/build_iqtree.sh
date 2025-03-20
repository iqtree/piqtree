echo "BOOST_ROOT=$BOOST_ROOT"
echo "Boost_INCLUDE_DIR=$Boost_INCLUDE_DIR"
echo "Boost_LIBRARY_DIRS=$Boost_LIBRARY_DIRS"
echo "CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH"

cd iqtree2
rm -rf build
mkdir build && cd build

if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Building for macOS."
    echo $LDFLAGS    
    echo $CPPFLAGS
    echo $CXXFLAGS    
    cmake -DBUILD_LIB=ON -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ ..
    gmake -j
elif [[ "$OSTYPE" == "msys"* || "$OSTYPE" == "cygwin"* ]]; then
    echo "Building for Windows."
    cmake -G "MinGW Makefiles" \
        -DCMAKE_C_COMPILER=$CMAKE_C_COMPILER \
        -DCMAKE_CXX_COMPILER=$CMAKE_CXX_COMPILER \
        -DCMAKE_MAKE_PROGRAM=$CMAKE_MAKE_PROGRAM \
        -DBUILD_LIB=ON \
        ..
    make -j
else
    echo "Building for Linux."
    cmake -DBUILD_LIB=ON ..
    make -j
fi

cd ../..

if [[ "$OSTYPE" == "darwin"* || "$OSTYPE" == "linux"* ]]; then
    mv iqtree2/build/libiqtree2.a src/piqtree/_libiqtree/
elif [[ "$OSTYPE" == "msys"* || "$OSTYPE" == "cygwin"* ]]; then
    mv iqtree2/build/iqtree2.lib src/piqtree/_libiqtree/
    mv iqtree2/build/iqtree2.dll src/piqtree/_libiqtree/
fi
