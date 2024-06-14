# Compile CarpetX (Frontier)

* Download CarpetX and SpacetimeX

    ```bash
    curl -kLO https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
    chmod a+x GetComponents
    ./GetComponents --root Cactus --parallel --no-shallow https://raw.githubusercontent.com/lwJi/SpacetimeX/main/Docs/thornlist/spacetimex.th
    ```

## The Short Way


### cce-17.0.0

* Load modules

```bash
source repos/SpacetimeX/Docs/compile-notes/frontier/Load-Module-CarpetX-cce17.sh
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
source repos/SpacetimeX/Docs/compile-notes/frontier/Load-Module-CarpetX-cce15.sh
```

* Install AsterX

```bash
cd Cactus
gmake AsterX options=simfactory/mdb/optionlists/frontier.cfg
cp repos/SpacetimeX/Docs/thornlist/asterx-frontier.th configs/AsterX/ThornList
gmake -j24 AsterX
```



## The Long Way

### Compile AMReX

* Clone `amrex`

```bash
git clone https://github.com/AMReX-Codes/amrex.git
```

* Install `amrex` to `$HOME/local/amrex-24.06`

```bash
cd amrex
mkdir build
cd build

source repos/SpacetimeX/Docs/compile-notes/frontier/amrex/Load-Module-AMReX.sh
source repos/SpacetimeX/Docs/compile-notes/frontier/amrex/Export-AMReX.sh
source repos/SpacetimeX/Docs/compile-notes/frontier/amrex/Compile-AMReX.sh

make -j24 install
```
