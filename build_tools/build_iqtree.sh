cd iqtree3
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

    if [[ -n "$BOOST_ROOT" ]]; then
        export Boost_INCLUDE_DIR="${BOOST_ROOT}"
        export Boost_LIBRARY_DIRS="${BOOST_ROOT}"
    fi

    cmake -G "MinGW Makefiles" \
        -DCMAKE_C_COMPILER=clang \
        -DCMAKE_CXX_COMPILER=clang++ \
        -DCMAKE_C_FLAGS=--target=x86_64-pc-windows-gnu \
        -DCMAKE_CXX_FLAGS=--target=x86_64-pc-windows-gnu \
        -DCMAKE_MAKE_PROGRAM=make \
        -DBoost_INCLUDE_DIR=$Boost_INCLUDE_DIR \
        -DBoost_LIBRARY_DIRS=$Boost_LIBRARY_DIRS \
        -DIQTREE_FLAGS="cpp14" \
        -DBUILD_LIB=ON \
        ..
    make -j
else
    echo "Building for linux."
    cmake -DBUILD_LIB=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..
    make -j
fi

cd ../..

if [[ "$OSTYPE" == "darwin"* || "$OSTYPE" == "linux"* ]]; then
    mv iqtree3/build/libiqtree.a src/piqtree/_libiqtree/
elif [[ "$OSTYPE" == "msys"* || "$OSTYPE" == "cygwin"* ]]; then
    mv iqtree3/build/iqtree.lib src/piqtree/_libiqtree/
    mv iqtree3/build/iqtree.dll src/piqtree/_libiqtree/
fi
