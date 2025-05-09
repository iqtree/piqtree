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
else
    echo "Building for linux."
    cmake -DBUILD_LIB=ON -DCMAKE_POLICY_VERSION_MINIMUM=3.5 ..
    make -j
fi

cd ../..
mv iqtree2/build/libiqtree2.a src/piqtree/_libiqtree/