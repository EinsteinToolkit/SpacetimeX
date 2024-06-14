# Compile on Frontier

## Compile CarpetX

### cce-17.0.0

* Load modules

```bash
source Load-Module-CarpetX-cce17.sh
```

* Install AsterX

```bash
cd Cactus
gmake AsterX options=repos/SpacetimeX/Docs/compile-notes/frontier/configs/frontier-cce17.cfg
cp repos/SpacetimeX/Docs/thornlist/asterx-frontier.th configs/AsterX/ThornList
gmake -j24 AsterX
```




### cce-15.0.0

* Load modules

```bash
source Load-Module-CarpetX-cce15.sh
```

* Install AsterX

```bash
cd Cactus
gmake AsterX options=simfactory/mdb/optionlists/frontier.cfg
cp repos/SpacetimeX/Docs/thornlist/asterx-frontier.th configs/AsterX/ThornList
gmake -j24 AsterX
```




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
