# \SCRIPT\-------------------------------------------------------------------------
#
#  CONTENTS      : download, build and install tools
#
#  DESCRIPTION   : none
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------
# Copyright (c) 2018-2020, Pay Giesselmann, Max Planck Institute for Molecular Genetics
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Written by Pay Giesselmann, modified by Gordon Rix
# ---------------------------------------------------------------------------------
# main build rules
import os, sys, site

rule default:
    shell : ""

rule all:
    input:
        "bin/minimap2",
        "bin/samtools",
        "bin/NGmerge",
        "lib/python3.8/C3POa.py",
        "lib/python3.8/site-packages/conk/conk.cpython-38-x86_64-linux-gnu.so",
        "lib/python3.8/site-packages/medaka/maple_smolecule.py"

rule all_but_guppy:
    input:
        "bin/minimap2",
        "bin/samtools",
        "bin/NGmerge",
        "lib/python3.8/C3POa.py",
        "lib/python3.8/site-packages/conk/conk.cpython-38-x86_64-linux-gnu.so",
        "lib/python3.8/site-packages/medaka/maple_smolecule.py"

# helper functions
def find_go():
    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, 'go')
        if os.path.isfile(exe_file):
            return exe_file
    return None

config['python'] = sys.executable

# defaults
if not 'threads_build' in config:
    config['threads_build'] = 1

if not 'build_generator' in config:
    config['build_generator'] = '"Unix Makefiles"'

if not 'flappie_src' in config:
    config['flappie_src'] = False

rule golang:
    output:
        go = "src/go/bin/go"
    shell:
        """
        mkdir -p src && cd src
        wget -q -nc https://dl.google.com/go/go1.13.4.linux-amd64.tar.gz
        tar -xzf go1.13.4.linux-amd64.tar.gz
        """

rule squashfs:
    output:
        mksquashfs = "bin/mksquashfs",
        unsquashfs = "bin/unsquashfs"
    shell:
        """
        install_prefix=`pwd`
        mkdir -p src && cd src
        if [ ! -d squashfs-tools ]; then
            git clone https://github.com/plougher/squashfs-tools --branch 4.4 --depth=1 && cd squashfs-tools
        else
            cd squashfs-tools && git fetch --all --tags --prune && git checkout tags/4.4
        fi
        cd squashfs-tools
        make
        cp *squashfs $install_prefix/bin/
        """

rule singularity:
    input:
        go = rules.golang.output.go
    output:
        "bin/singularity"
    shell:
        """
        install_prefix=`pwd`
        mkdir -p src/gocode
        export GOPATH=$(pwd)/src/gocode
        export GOROOT=$(pwd)/src/go
        export PATH=$(pwd)/$(dirname {input.go}):$PATH
        cd src
        if [ ! -d singularity ]; then
            git clone https://github.com/sylabs/singularity.git --branch v3.3.0 --depth=1 && cd singularity
        else
            cd singularity && git fetch --all --tags --prune && git checkout tags/v3.3.0
        fi
        ./mconfig --without-suid --prefix=$install_prefix --localstatedir=$install_prefix/
        make -C ./builddir
        make -C ./builddir install
        """

rule vbz_compression:
    output:
        lib = "lib/libvbz_hdf_plugin.so.1.0.1"
    shell:
        """
        install_prefix=`pwd`
        mkdir -p lib
        mkdir -p src && cd src
        if [ ! -d singularity ]; then
            git clone https://github.com/nanoporetech/vbz_compression --recursive --branch v1.0.1 --depth=1 && cd vbz_compression
        else
            cd vbz_compression && git fetch --all --tags --prune && git checkout tags/v1.0.1
        fi
        rm -rf build
        mkdir build && cd build
        cmake -D CMAKE_BUILD_TYPE=Release -D ENABLE_CONAN=OFF -D ENABLE_PERF_TESTING=OFF -D ENABLE_PYTHON=OFF ..
        make
        cp bin/libvbz_hdf_plugin.so ${{install_prefix}}/{output.lib}
        """

# detailed build rules
rule htslib:
    output:
        src = directory("src/htslib")
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d htslib ]; then
            git clone https://github.com/samtools/htslib --branch 1.9 --depth=1 && cd htslib
        else
            cd htslib && git fetch --all --tags --prune && git checkout tags/1.9
        fi
        autoheader && autoconf && ./configure && make -j{threads}
        """

rule samtools:
    input:
        rules.htslib.output.src
    output:
        bin = "bin/samtools"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d samtools ]; then
            git clone https://github.com/samtools/samtools --branch 1.9 --depth=1 && cd samtools
        else
            cd samtools && git fetch --all --tags --prune && git checkout tags/1.9
        fi
        autoheader --warning=none && autoconf -Wno-syntax && ./configure && make -j{threads}
        cp samtools ../../{output.bin}
        """

rule minimap2:
    output:
        bin = "bin/minimap2"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d minimap2 ]; then
            git clone https://github.com/lh3/minimap2 --branch v2.14 --depth=1 && cd minimap2
        else
            cd minimap2 && git fetch --all --tags --prune && git checkout tags/v2.14
        fi
        make clean && make -j{threads}
        cp minimap2 ../../{output.bin}
        """


rule NGmerge:
    output:
        bin = "bin/NGmerge"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d NGmerge ]; then
            git clone https://github.com/harvardinformatics/NGmerge --branch v0.3 --depth=1 && cd NGmerge
        else
            cd NGmerge && git fetch --all --tags --prune && git checkout tags/v0.3
        fi
        make clean && make -j{threads}
        cp NGmerge ../../{output.bin}
        """

rule C3POa:
    output:
        C3POa = "lib/python3.8/C3POa.py",
        conk = "lib/python3.8/site-packages/conk/conk.cpython-38-x86_64-linux-gnu.so"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ -d C3POa ]; then
            rm -r -f C3POa
        fi
        git clone -b peakFinderCustomSettings https://github.com/christopher-vollmers/C3POa.git

        if [ -d conk ]; then
            rm -r -f conk
        fi
        git clone https://github.com/rvolden/conk && cd conk
        python setup.py sdist bdist_wheel
        python -m pip install dist/conk*whl --force-reinstall
        cd ../..
        
        if [ -d lib/python3.8/bin ]; then
            mv src/C3POa/bin/* lib/python3.8/bin/
        else
            mv src/C3POa/bin lib/python3.8
        fi
        mv src/C3POa/C3POa.py src/C3POa/C3POa_postprocessing.py lib/python3.8
        rm -r -f src/C3POa
        """

rule maple_medaka:
    output:
        ms = "lib/python3.8/site-packages/medaka/maple_smolecule.py"
    shell:
        """
        mkdir -p src && cd src
        if [ -d maple ]; then
            rm -r -f maple
        fi
        git clone https://github.com/gordonrix/maple.git && cd ..
        mv src/maple/rules/utils/maple_smolecule.py src/maple/rules/utils/medaka.py lib/python3.8/site-packages/medaka
        rm -r -f src/maple
        """


rule UCSCtools:
    output:
        "bin/bedGraphToBigWig"
    shell:
        """
        wget -q -O {output} https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
        chmod 755 {output}
        """

rule bedtools:
    output:
        bin = "bin/bedtools"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d bedtools2 ]; then
            git clone https://github.com/arq5x/bedtools2 --branch v2.27.1 --depth=1 && cd bedtools2
        else
            cd bedtools2 && git fetch --all --tags --prune && git checkout tags/v2.27.1
        fi
        make clean && make
        cp bin/bedtools ../../{output.bin}
        """


rule ngmlr:
    output:
        bin = "bin/ngmlr"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d ngmlr ]; then
            git clone https://github.com/philres/ngmlr --branch v0.2.7 --depth=1 && cd ngmlr
        else
            cd ngmlr && git fetch --all --tags --prune && git checkout tags/v0.2.7
        fi
        mkdir -p build && cd build && rm -rf * && cmake -DCMAKE_BUILD_TYPE=Release -G {config[build_generator]} .. && cmake --build . --config Release -- -j {threads}
        cp ../bin/*/ngmlr ../../../bin
        """

rule nanopolish:
    input:
        rules.vbz_compression.output.lib
    output:
        bin = "bin/nanopolish"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d nanopolish ]; then
            #git clone --recursive https://github.com/jts/nanopolish --branch v0.11.0 --depth=1 && cd nanopolish
            git clone --recursive https://github.com/jts/nanopolish && cd nanopolish && git reset --hard 46a029c
        else
            #cd nanopolish && git fetch --all --tags --prune && git checkout tags/v0.11.0
            cd nanopolish && git fetch --all --tags --prune && git checkout 46a029c
        fi
        make clean
        set +e
        make -j{threads} 2>&1 | tee log.txt | grep 'g++\|gcc'
        exitcode=$?
        if [ $exitcode -ne 0 ]; then tail -n 100 log.txt; exit 1; fi
        cp nanopolish ../../bin/
        """

rule sniffles:
    output:
        bin = "bin/sniffles"
    threads: config['threads_build']
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d Sniffles ]; then
            git clone https://github.com/fritzsedlazeck/Sniffles --branch v1.0.10 --depth=1 && cd Sniffles
        else
            cd Sniffles && git fetch --all --tags --prune && git checkout tags/v1.0.10
        fi
        mkdir -p build && cd build && rm -rf * && cmake -DCMAKE_BUILD_TYPE=Release -G{config[build_generator]} .. && cmake --build . --config Release -- -j {threads}
        cp ../bin/*/sniffles ../../../bin
        """

rule svim:
    output:
        bin = "bin/svim"
    params:
        prefix = lambda wildcards : sys.prefix
    shell:
        """
        prefix=`pwd`
        mkdir -p src && cd src
        if [ ! -d svim ]; then
            git clone https://github.com/eldariont/svim.git --branch v1.4.0 --depth=1 && cd svim
        else
            cd svim && git fetch --all --tags --prune && git checkout tags/v1.4.0
        fi
        {config[python]} -m pip install .
        find {params.prefix}/bin {params.prefix}/local/bin -type f -name svim -exec cp {{}} ../../{output.bin} \; || true
        """

rule deepbinner:
    input:
        rules.vbz_compression.output.lib
    output:
        bin = 'bin/deepbinner'
    params:
        prefix = lambda wildcards : sys.prefix
    shell:
        """
        mkdir -p src && cd src
        if [ ! -d Deepbinner ]; then
            git clone https://github.com/rrwick/Deepbinner --branch v0.2.0 --depth=1 && cd Deepbinner
        else
            cd Deepbinner && git fetch --all --tags --prune && git checkout tags/v0.2.0
        fi
        {config[python]} -m pip install "tensorflow==1.15" "keras==2.2.5"
        {config[python]} -m pip install .
        find {params.prefix}/bin {params.prefix}/local/bin -type f -name deepbinner -exec cp {{}} ../../{output.bin} \; || true
        """

rule gitlfs:
    input:
        go = lambda wildcards : find_go() if find_go() is not None else rules.golang.output.go
    output:
        bin = "bin/git-lfs"
    shell:
        """
        mkdir -p src/gocode
        export GOPATH=$(pwd)/src/gocode
        {input.go} get github.com/github/git-lfs
        cp src/gocode/bin/git-lfs bin/
        """

rule OpenBLAS:
    output:
        src = directory("src/OpenBLAS")
    threads: config['threads_build']
    shell:
        """
        install_prefix=`pwd`
        mkdir -p src && cd src
        if [ ! -d OpenBLAS ]; then
            git clone https://github.com/xianyi/OpenBLAS --branch v0.3.4 --depth=1 && cd OpenBLAS
        else
            cd OpenBLAS && git fetch --all --tags --prune && git checkout tags/v0.3.4
        fi
        make clean && make NO_LAPACK=1 NOFORTRAN=1 -j{threads}
        make install PREFIX=$install_prefix
        """

rule hdf5:
    output:
        src = directory("src/hdf5")
    threads: config['threads_build']
    shell:
        """
        install_prefix=`pwd`
        mkdir -p src && cd src
        if [ ! -d hdf5 ]; then
            git clone https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git --branch hdf5-1_8_20 --depth=1 && cd hdf5
        else
            cd hdf5 && git fetch --all --tags --prune && git checkout tags/hdf5-1_8_20
        fi
        mkdir -p build && cd build && rm -rf * && cmake -DCMAKE_BUILD_TYPE=Release -DHDF5_ENABLE_Z_LIB_SUPPORT=ON -DHDF5_BUILD_TOOLS=OFF -DBUILD_TESTING=OFF -DHDF5_BUILD_EXAMPLES=OFF -DCMAKE_INSTALL_PREFIX=$install_prefix -G{config[build_generator]} ../
        cmake --build . --config Release -- -j {threads}
        cmake --build . --config Release --target install
        """

