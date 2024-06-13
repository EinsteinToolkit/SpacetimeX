# Compile on Frontier

## Compile AMReX

```bash
cd amrex
mkdir build
cd build

source Load-Module-AMReX.sh
source Export-AMReX.sh
source Compile-AMReX.sh

make -j8 install
make test_install
```
