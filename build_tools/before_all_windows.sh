# Install dependencies using choco
choco install -y llvm --version=14.0.6 --allow-downgrade
choco install -y eigen 

# Build IQ-TREE
bash build_tools/build_iqtree.sh