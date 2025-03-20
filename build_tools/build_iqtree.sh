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
        -DCMAKE_C_COMPILER=clang \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DCMAKE_C_FLAGS=--target=x86_64-pc-windows-gnu \
        -DCMAKE_CXX_FLAGS=--target=x86_64-pc-windows-gnu \
        -DCMAKE_MAKE_PROGRAM=mingw32-make \
        -DBoost_INCLUDE_DIR=$Boost_INCLUDE_DIR \
        -DBoost_LIBRARY_DIRS=$Boost_LIBRARY_DIRS \
        -DIQTREE_FLAGS="cpp14" \
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
