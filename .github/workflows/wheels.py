# Workflow to build and test wheels
name: Wheel builder

on:
  schedule:
    # Nightly build at 3:42 A.M.
    - cron: "42 3 */1 * *"
  push:
    branches:
      # - main
      # Feature branch
      - feature-cibuildwheel
      # Release branches
      - "[0-9]+.[0-9]+.X"
  pull_request:
    branches:
      - main
      - "[0-9]+.[0-9]+.X"
  # Manual run
  workflow_dispatch:

concurrency:
  group: ${{ github.workflow }}-${{ github.head_ref || github.run_id }}
  cancel-in-progress: true

jobs:
  # Check whether to build the wheels and the source tarball
  check_build_trigger:
    name: Check build trigger
    runs-on: ubuntu-latest
    if: github.repository == 'UG4/ugcore'
    outputs:
      build: ${{ steps.check_build_trigger.outputs.build }}

    steps:
      - name: Checkout UG4
        uses: actions/checkout@v3
      #  with:
      #    ref: ${{ github.event.pull_request.head.sha }}

      #- id: check_build_trigger
      #  name: Check build trigger
      #  run: bash build_tools/github/check_build_trigger.sh

  # Build the wheels for Linux, Windows and macOS for Python 3.9 and newer
  build_wheels:
    name: Build wheel for cp${{ matrix.python }}-${{ matrix.platform_id }}-${{ matrix.manylinux_image }}
    runs-on: ${{ matrix.os }}
    needs: check_build_trigger
    if: needs.check_build_trigger.outputs.build

    strategy:
      # Ensure that a wheel builder finishes even if another fails
      fail-fast: false
      matrix:
        include:
          # Window 64 bit
          # Note: windows-2019 is needed for older Python versions:
          # https://github.com/scikit-learn/scikit-learn/issues/22530
          - os: windows-latest
            python: 39
            platform_id: win_amd64
          - os: windows-latest
            python: 310
            platform_id: win_amd64
          - os: windows-latest
            python: 311
            platform_id: win_amd64
          - os: windows-latest
            python: 312
            platform_id: win_amd64

          # Linux 64 bit manylinux2014
          - os: ubuntu-latest
            python: 39
            platform_id: manylinux_x86_64
            manylinux_image: manylinux2014

          # NumPy on Python 3.10 only supports 64bit and is only available with manylinux2014
          - os: ubuntu-latest
            python: 310
            platform_id: manylinux_x86_64
            manylinux_image: manylinux2014

          - os: ubuntu-latest
            python: 311
            platform_id: manylinux_x86_64
            manylinux_image: manylinux2014
          - os: ubuntu-latest
            python: 312
            platform_id: manylinux_x86_64
            manylinux_image: manylinux2014

          # MacOS x86_64
          - os: macos-latest
            python: 39
            platform_id: macosx_x86_64
          - os: macos-latest
            python: 310
            platform_id: macosx_x86_64
          - os: macos-latest
            python: 311
            platform_id: macosx_x86_64
          - os: macos-latest
            python: 312
            platform_id: macosx_x86_64

          # MacOS arm64
          - os: macos-14
            python: 39
            platform_id: macosx_arm64
          - os: macos-14
            python: 310
            platform_id: macosx_arm64
          - os: macos-14
            python: 311
            platform_id: macosx_arm64
          - os: macos-14
            python: 312
            platform_id: macosx_arm64

    steps:
      - name: Checkout scikit-learn
        uses: actions/checkout@v3

      - name: Setup Python
        uses: actions/setup-python@v5
        with:
          python-version: "3.11" # update once build dependencies are available

      #  run: bash build_tools/wheels/build_wheels.sh
     - name: Install cibuildwheel
        run: python -m pip install cibuildwheel==2.17.0

      - name: Build wheels
        run: python -m cibuildwheel --output-dir wheelhouse
        # to supply options, put them in 'env', like:
        # env:
         #    CIBW_PRERELEASE_PYTHONS: ${{ matrix.prerelease }}
      #    CIBW_ENVIRONMENT: SKLEARN_SKIP_NETWORK_TESTS=1
      #      SKLEARN_BUILD_PARALLEL=3
      #    CIBW_BUILD: cp${{ matrix.python }}-${{ matrix.platform_id }}
      #    CIBW_ARCHS: all
      #    CIBW_MANYLINUX_X86_64_IMAGE: ${{ matrix.manylinux_image }}
      #    CIBW_MANYLINUX_I686_IMAGE: ${{ matrix.manylinux_image }}
      #    # Needed on Windows CI to compile with Visual Studio compiler
      #    # otherwise Meson detects a MINGW64 platform and use MINGW64
      #    # toolchain
      #    CIBW_CONFIG_SETTINGS_WINDOWS: "setup-args=--vsenv"
      #    CIBW_REPAIR_WHEEL_COMMAND_WINDOWS: bash build_tools/github/repair_windows_wheels.sh {wheel} {dest_dir}
      #    CIBW_BEFORE_TEST_WINDOWS: bash build_tools/github/build_minimal_windows_image.sh ${{ matrix.python }}
      #    CIBW_TEST_REQUIRES: pytest pandas
      #    CIBW_TEST_COMMAND: bash {project}/build_tools/wheels/test_wheels.sh
      #    CIBW_TEST_COMMAND_WINDOWS: bash {project}/build_tools/github/test_windows_wheels.sh ${{ matrix.python }}
      #    CIBW_BUILD_VERBOSITY: 1

      - name: Store artifacts
        uses: actions/upload-artifact@v3
        with:
          path: wheelhouse/*.whl
