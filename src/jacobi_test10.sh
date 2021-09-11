rm myspkmeanssp.cpython-37m-x86_64-linux-gnu.so
python3 setup.py build_ext --inplace
python3 spkmeans.py 0 jacobi ../tests/matrix10.csv