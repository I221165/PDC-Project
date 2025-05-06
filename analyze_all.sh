#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <k>"
  exit 1
fi
K=$1

# dataset paths, relative to each build-*/ dir
declare -A EDGEFILES=(
  [10k]="../data/out2.edgelist"
  [100k]="../data/sample10k.edgelist"
  [full]="../data/higgs/out.edgelist"
)
declare -A INTFILES=(
  [10k]="../data/interests2.csv"
  [100k]="../data/interests10k.csv"
  [full]="../data/higgs/interests.csv"
)

SIZES=(10k 100k full)
THREADS=(1 4)   # ensure this is defined before both loops

########################################
# 1) Build with -pg for gprof
########################################
echo "=== Building for gprof ==="
rm -rf build-gprof
mkdir build-gprof && cd build-gprof
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_FLAGS="-O3 -pg" \
      -DCMAKE_EXE_LINKER_FLAGS="-pg" \
      ..
make -j
cd ..

########################################
# 2) gprof runs over sizes × threads
########################################
for size in "${SIZES[@]}"; do
  for th in "${THREADS[@]}"; do
    echo ">>> gprof run on ${size} with OMP_NUM_THREADS=${th}"
    pushd build-gprof >/dev/null
      export OMP_NUM_THREADS=$th
      mpirun -np 1 ./psaiim "${EDGEFILES[$size]}" "${INTFILES[$size]}" "$K"
      mv gmon.out gmon_${size}_Th${th}.out
      gprof ./psaiim gmon_${size}_Th${th}.out > profile_${size}_Th${th}.txt
      echo "-> profile_${size}_Th${th}.txt"
    popd >/dev/null
  done
done

########################################
# 3) Build with -g for Valgrind
########################################
echo "=== Building for Valgrind ==="
rm -rf build-valgrind
mkdir build-valgrind && cd build-valgrind
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j
cd ..

########################################
# 4) Valgrind runs over sizes × threads
########################################
for size in "${SIZES[@]}"; do
  for th in "${THREADS[@]}"; do
    echo ">>> Valgrind run on ${size} with OMP_NUM_THREADS=${th}"
    pushd build-valgrind >/dev/null
      export OMP_NUM_THREADS=$th
      mpirun -np 1 valgrind \
        --leak-check=full \
        --show-leak-kinds=all \
        --track-origins=yes \
        --log-file=valgrind_${size}_Th${th}.log \
        ./psaiim "${EDGEFILES[$size]}" "${INTFILES[$size]}" "$K"
      echo "-> valgrind_${size}_Th${th}.log"
    popd >/dev/null
  done
done
