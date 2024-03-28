# Install CarpetX with Spack (Frontera)

* Use interactive session

    - Compile CPU version: `idev -m 120`

    - Compile GPU version: `idev -p rtx-dev -m 120`

* Download CarpetX and SpacetimeX

    ```bash
    curl -kLO https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
    chmod a+x GetComponents
    ./GetComponents --root Cactus --parallel --no-shallow https://raw.githubusercontent.com/lwJi/SpacetimeX/main/Docs/thornlist/spacetimex.th
    ```


## The Short Way

### Intel-Oneapi version (`oneapi@2023.1.0`)

* Load `intel/23.1.0`

    - `module load intel/23.1.0`

* Install SpacetimeX

    ```bash
    cd Cactus
    gmake SpacetimeX-oneapi options=repos/SpacetimeX/Docs/compile-notes/frontera/configs/config-frontera-oneapi.cfg
    cp repos/SpacetimeX/Docs/thornlist/spacetimex.th configs/SpacetimeX-oneapi/ThornList
    gmake -j24 SpacetimeX-oneapi
    ```

### CUDA version (`cuda@11.8.0`)

* Load `intel/19.1.1`

    - `module reset` or `module load intel/19.1.1`

* Install SpacetimeX

    ```bash
    . /work2/08708/liwei/frontera/SpackSource/spack/share/spack/setup-env.sh
    spack load gcc@11.2.0
    spack load cuda@11.8.0

    cd Cactus
    gmake SpacetimeX-cuda options=repos/SpacetimeX/Docs/compile-notes/frontera/configs/config-frontera-gcc-cuda-impi.cfg
    cp repos/SpacetimeX/Docs/thornlist/spacetimex.th configs/SpacetimeX-cuda/ThornList
    gmake -j16 SpacetimeX-cuda
    ```


## The Long Way

### Intel-Oneapi version (`oneapi@2023.1.0`)

* Load `intel/23.1.0`

    - `module load intel/23.1.0`

* Download spack

    - `git clone -c feature.manyFiles=true https://github.com/spack/spack.git`

    - use `develop` branch: (maybe `git checkout 141c7de5`)

    - `. share/spack/setup-env.sh`
    
    - `spack compiler find`

* Create a dir where you want put `view` in (say `/work2/.../username/frontera/SpackView/oneapi`)

* Replace the last line of `oneapi-23.1.0/spack.yaml` with your own dir (say `/work2/.../username/frontera/SpackView/oneapi/view`)

* Replace the dir `/work2/08708/liwei/frontera/SpackView/oneapi/view` (with say `/work2/.../username/frontera/SpackView/oneapi/view`) in `configs/config-frontera-oneapi.cfg`

* Install other required packages

    ```bash
    env TMPDIR=$WORK/tmp spack --env-dir ./oneapi-23.1.0 compiler find
    env TMPDIR=$WORK/tmp spack --env-dir ./oneapi-23.1.0 concretize --force
    env TMPDIR=$WORK/tmp spack --env-dir ./oneapi-23.1.0 install --fail-fast
    ```

* Install SpacetimeX

    ```bash
    cd Cactus
    gmake SpacetimeX-oneapi options=repos/SpacetimeX/Docs/compile-notes/frontera/configs/config-frontera-oneapi.cfg
    cp repos/SpacetimeX/Docs/thornlist/spacetimex.th configs/SpacetimeX-oneapi/ThornList
    gmake -j24 SpacetimeX-oneapi
    ```


### Cuda version (`cuda@11.8.0` with `intel@23.1.0`)

* Load `intel/19.1.1`

    - `module reset` or `module load intel/19.1.1`

* Download spack

    - `git clone -c feature.manyFiles=true https://github.com/spack/spack.git`

    - `git checkout relesases/v0.21`

    - `. share/spack/setup-env.sh`

    - `spack compiler find`

* Install gcc@11.2.0

```bash
spack install gcc@11.2.0 %gcc@4.8.5
spack compiler add ...  # ... is the last line of previous command
```

* Create a dir where you want put `view` in (say `/work2/.../username/frontera/SpackView/gcc11.2.0-cuda11.8.0-impi19.0.9`)

* Replace the last line of `gcc11.2.0-cuda11.8.0-impi19.0.9/spack_yaml` with your dir (say `/work2/.../username/frontera/SpackView/gcc11.2.0-cuda11.8.0-impi19.0.9/view`)

* Replace the dir `/work2/08708/liwei/frontera/SpackView/gcc11.2.0-cuda11.8.0-impi19.0.9/view` (with say `/work2/.../username/frontera/SpackView/gcc11.2.0-cuda11.8.0-impi19.0.9/view`)
in `config-frontera-gcc11.2.0-cuda11.8.0-impi19.0.9.cfg`

* Install other required packages

```bash
env TMPDIR=$WORK/tmp spack --env-dir ./gcc11.2.0-cuda11.8.0-impi19.0.9 compiler find view-cuda-compilers
env TMPDIR=$WORK/tmp spack --env-dir ./gcc11.2.0-cuda11.8.0-impi19.0.9 concretize --force
env TMPDIR=$WORK/tmp spack --env-dir ./gcc11.2.0-cuda11.8.0-impi19.0.9 install --fail-fast
```

* Install SpacetimeX

```bash
spack load gcc@11.2.0
spack load cuda@11.8.0
cd Cactus
gmake SpacetimeX-cuda options=config-frontera-gcc-cuda-impi.cfg
cp repos/SpacetimeX/Docs/thornlist/spacetimex.th configs/SpacetimeX-cuda/ThornList
gmake -j16 SpacetimeX-cuda
```
