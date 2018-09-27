# Note

**Code review of *[IFOPT](https://github.com/ethz-adrl/ifopt)*, thanks for their work!**

# Build & Test

```shell
# prepare dependencies
sudo apt-get install cmake libeigen3-dev coinor-libipopt-dev

# build
mkdir build && cd build
cmake ..
make -j2

# run test
../bin/ex_test_ipopt

```
