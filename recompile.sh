cd $ISSM_DIR
autoreconf -ivf 
./configure.sh
make -j 8 install
