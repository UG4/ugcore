name: Build

on: [push, pull_request]

jobs:
  build_wheels:
    name: Build wheels on ${{ matrix.os }}
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, windows-latest, macos-13, macos-14]

    steps:
      - uses: actions/checkout@v4

      # Used to host cibuildwheel
      - uses: actions/setup-python@v3

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

      - uses: actions/upload-artifact@v4
        with:
          name: cibw-wheels-${{ matrix.os }}-${{ strategy.job-index }}
          path: ./wheelhouse/*.whl
