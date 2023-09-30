# Install CarpetX with Spack (Frontera)

Use interactive session

* Compile CPU version: `idev -m 120`

* Compile GPU version: `idev -p rtx-dev -m 120`

Download spack

* `git clone -c feature.manyFiles=true https://github.com/spack/spack.git`

* `git checkout relesases/v0.20`

## Download CarpetX and SpacetimeX

```
curl -kLO https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
chmod a+x GetComponents
./GetComponents --root Cactus --parallel --no-shallow https://raw.githubusercontent.com/lwJi/SpacetimeX/main/scripts/spacetimex.th
```

```
cd Cactus/repos
git clone https://github.com/lwJi/SpacetimeX.git
cd ../arrangements
ln -s ../repos/SpacetimeX
```

## Intel-Oneapi version (`oneapi@2023.1.0`)

Load `intel/23.1.0`

* `module load intel/23.1.0`

Setup spack

* `. share/spack/setup-env.sh`

* `spack compiler find`

Create a dir where you want put `view` in (say `/work2/.../username/frontera/SpackView/oneapi`)

* replace the last line of `oneapi-23.1.0/spack.yaml` with your dir (say `/work2/.../username/frontera/SpackView/oneapi/view`)

* replace the dir `/work2/08708/liwei/frontera/SpackView/oneapi/view` (with say `/work2/.../username/frontera/SpackView/oneapi/view`)
in `config_frontera_oneapi.cfg`

Install other required packages

* `env TMPDIR=$WORK/tmp spack --env-dir ./oneapi-23.1.0 compiler find`

* `env TMPDIR=$WORK/tmp spack --env-dir ./oneapi-23.1.0 concretize --force`

* `env TMPDIR=$WORK/tmp spack --env-dir ./oneapi-23.1.0 install --fail-fast`

Install CarpetX

* `cd Cactus`

* `gmake CarpetX-oneapi options=config_frontera_oneapi.cfg`

* `cp repos/AsterX/scripts/asterx.th configs/CarpetX-oneapi/ThornList`

* `gmake -j16 CarpetX-oneapi`


