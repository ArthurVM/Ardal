python -m pip uninstall -y ardal || true
rm -rf build/ dist/ *.egg-info

export CXXFLAGS="-O0 -g -fno-omit-frame-pointer -D_GLIBCXX_ASSERTIONS"
export CFLAGS="-O0 -g -fno-omit-frame-pointer -D_GLIBCXX_ASSERTIONS"
export LDFLAGS="-g"
export OMP_NUM_THREADS=1

python -m pip install -v -e .