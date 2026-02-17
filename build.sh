#!/bin/bash

set -euo pipefail

help()
{
   echo "Basic script to build mercator"
   echo
   echo "Syntax: ./build.sh [-b|-s]"
   echo "options:"
   echo "h                 Print this help."
   echo "b [BUILD_TYPE]    Build type: {Release, Debug}."
   echo "s [SYSTEM]        OS type: {win, linux, mac}."
   echo
}

while getopts ":hb:s:" option; do
   case $option in
      h)
        help
        exit;;
      b)
        build_type=$OPTARG;;
      s)
        os_type=$OPTARG;;
     \?)
        echo "Error: Invalid option"
        exit;;
   esac
done

if [[ -z "${build_type:-}" ]]; then
   build_type="Release"
fi

if [[ -z "${os_type:-}" ]]; then
   case "$(uname -s)" in
      Darwin)
         os_type="mac"
         ;;
      *)
         os_type="linux"
         ;;
   esac
fi

mkdir -p build && cd build/ # go to build folder
#export OMP_NUM_THREADS=$(nproc)

# Windows 10
if [[ $os_type == "win" ]]; then
    cmake .. -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=$build_type
    cmake --build . --config $build_type -j 8
fi

# Linux or Mac
if [[ "$os_type" == "linux" || "$os_type" == "mac" ]]; then
    cmake_args=(.. -G "Unix Makefiles" -DCMAKE_BUILD_TYPE="$build_type")

    if [[ "$os_type" == "mac" ]]; then
        if [[ -n "${SDKROOT:-}" && ! -d "${SDKROOT}" ]]; then
            echo "Warning: SDKROOT points to a missing path (${SDKROOT}); unsetting it for this build."
            unset SDKROOT
        fi

        macos_sdk_path=""
        if command -v xcrun >/dev/null 2>&1; then
            macos_sdk_path="$(xcrun --sdk macosx --show-sdk-path 2>/dev/null || true)"
        fi

        # Override stale cached sysroots from older Xcode/CLT installations.
        cmake_args+=("-DCMAKE_OSX_SYSROOT=${macos_sdk_path}")
    fi

    cmake "${cmake_args[@]}"
    cmake --build . -j 8
fi

if [[ -f bmercator ]]; then
    mv bmercator ../ # copy mercator to project directory
else
    echo "Error: build completed without producing ./build/bmercator."
    exit 1
fi
