#!/bin/bash
# Build SBNSoftware/sbn_ml_cafmaker (SPINE HDF5 -> CAF merging) against the
# sbnanaobj of the campaign's icaruscode release, in the SL7 container:
#
#   apptainer exec --cleanenv -B /cvmfs,/exp,/nashome,/pnfs/icarus \
#       /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest \
#       bash batch/ncdeltarad/build_ml_cafmaker.sh /exp/icarus/app/users/$USER/ml_cafmaker
#
# Arguments: <prefix> [icaruscode version, default v10_06_00_13] [qual, default e26:prof]
#
# Steps (each re-run skips what is already built):
#   1. sbnanaobj version = the one icaruscode <version> depends on (the CAFs
#      were written with it); sbnana = the newest sbnana built on that same
#      sbnanaobj (flatten_caf must read what the merge writes);
#   2. HDF5 1.10.5 with the C++ API from source (UPS hdf5 has no C++ API for
#      e26 -- sbn_ml_cafmaker README);
#   3. sbn_ml_cafmaker merge_sources_simulation;
#   4. <prefix>/setup_ml_cafmaker.sh: the runtime setup, read by
#      spine_caf.py --cafmaker-dir <prefix>.
# No `set -e`: the UPS setup scripts return non-zero internally.

PREFIX=${1:?usage: build_ml_cafmaker.sh <prefix> [icaruscode version] [qual]}
REL=${2:-v10_06_00_13}
QUAL=${3:-e26:prof}
die() { echo "build_ml_cafmaker: FAILED: $*" >&2; exit 1; }

mkdir -p "$PREFIX" && cd "$PREFIX" || die "cannot use $PREFIX"
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh || die "setup_icarus.sh"

OBJ=$(ups depend icaruscode "$REL" -q "$QUAL" 2>/dev/null | grep -o 'sbnanaobj v[0-9a-z_]*' | head -1 | cut -d' ' -f2)
[ -n "$OBJ" ] || die "no sbnanaobj in 'ups depend icaruscode $REL -q $QUAL'"
SBNANA=
for v in $(ups list -aK VERSION sbnana -q "$QUAL" 2>/dev/null | tr -d '"' | sort -rV | uniq); do
    if ups depend sbnana "$v" -q "$QUAL" 2>/dev/null | grep -q "sbnanaobj $OBJ\b"; then SBNANA=$v; break; fi
done
[ -n "$SBNANA" ] || die "no sbnana -q $QUAL depends on sbnanaobj $OBJ (ups list -aK+ sbnana)"
echo "icaruscode $REL -> sbnanaobj $OBJ -> sbnana $SBNANA"
setup sbnana "$SBNANA" -q "$QUAL" || die "setup sbnana $SBNANA -q $QUAL"
command -v cmake >/dev/null || setup cmake || die "no cmake"
[ -f "$SBNANAOBJ_INC/sbnanaobj/StandardRecord/SRInteractionDLP.h" ] \
    || die "sbnanaobj $OBJ has no SRInteractionDLP (the ML CAF branches)"

# 2. HDF5 C++
H=$PREFIX/hdf5
if [ ! -e "$H/install/lib/libhdf5_cpp.so" ]; then
    [ -d "$H/src" ] || git clone -q --depth 1 --branch hdf5-1_10_5 https://github.com/HDFGroup/hdf5.git "$H/src" \
        || die "git clone hdf5"
    cmake -S "$H/src" -B "$H/build" -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CXX_STANDARD=17 \
          -DCMAKE_INSTALL_PREFIX="$H/install" -DHDF5_BUILD_CPP_LIB=ON -DHDF5_BUILD_HL_LIB=ON \
          -DBUILD_TESTING=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_TOOLS=OFF >"$H/cmake.log" 2>&1 \
        || die "hdf5 cmake (see $H/cmake.log)"
    cmake --build "$H/build" -j 8 --target install >"$H/build.log" 2>&1 || die "hdf5 build (see $H/build.log)"
fi

# 3. sbn_ml_cafmaker
S=$PREFIX/sbn_ml_cafmaker
if [ -d "$S/.git" ]; then git -C "$S" pull -q || die "git pull sbn_ml_cafmaker"
else git clone -q https://github.com/SBNSoftware/sbn_ml_cafmaker.git "$S" || die "git clone sbn_ml_cafmaker"; fi
cmake -S "$S" -B "$S/build" -DHDF5_INSTALL="$H/install" >"$S/cmake.log" 2>&1 || die "cmake (see $S/cmake.log)"
cmake --build "$S/build" -j 8 --target merge_sources_simulation >"$S/build.log" 2>&1 \
    || { tail -30 "$S/build.log" >&2; die "build (see $S/build.log)"; }

# 4. runtime setup + smoke test
cat > "$PREFIX/setup_ml_cafmaker.sh" <<EOF
source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana $SBNANA -q $QUAL
export LD_LIBRARY_PATH=$H/install/lib:\$LD_LIBRARY_PATH
EOF
export LD_LIBRARY_PATH=$H/install/lib:$LD_LIBRARY_PATH
missing=$(ldd "$S/build/merge_sources_simulation" | grep 'not found')
[ -z "$missing" ] || die "unresolved libraries: $missing"
command -v flatten_caf >/dev/null || die "no flatten_caf in sbnana $SBNANA"
echo "built $S/build/merge_sources_simulation ($(git -C "$S" log -1 --format='%h %ci'))"
echo "sbnana $SBNANA / sbnanaobj $OBJ; runtime setup: $PREFIX/setup_ml_cafmaker.sh"
echo "use: spine_caf.py <workdir> --h5 '<h5 glob>' --cafmaker-dir $PREFIX"
