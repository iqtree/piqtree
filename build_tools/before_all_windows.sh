# Install dependencies using choco

echo "Boost_INCLUDE_DIR: " $Boost_INCLUDE_DIR
echo "Boost_LIBRARY_DIRS: " $Boost_LIBRARY_DIRS
echo "BOOST_ROOT: " $BOOST_ROOT

choco install -y llvm --version=14.0.6 --allow-downgrade
choco install -y eigen 

# Build IQ-TREE
bash build_tools/build_iqtree.sh