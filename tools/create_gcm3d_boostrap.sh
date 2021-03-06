#!/bin/bash

SCRIPT="`readlink -f $0`"
DIR="`pwd`"

get_attachment() {
    sed -e "1,/^# ----- attachment $1/d" "$SCRIPT" | sed -e '/^# ----- attachment /,$d'
}

get_build_script() {
    get_attachment 1
}

get_bootstrap_script() {
    get_attachment 2
}

cd /tmp
rm -rf gcm-3d-bundle && mkdir gcm-3d-bundle && cd gcm-3d-bundle

# Links to original sources (as for 08.07.2013)
#
# http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.gz
# ftp://ftp.gmplib.org/pub/gmp/gmp-5.1.2.tar.bz2
# http://www.multiprecision.org/mpc/download/mpc-1.0.1.tar.gz
# ftp://gcc.gnu.org/pub/gcc/infrastructure/cloog-0.18.0.tar.gz
# ftp://gcc.gnu.org/pub/gcc/infrastructure/isl-0.11.1.tar.bz2
# http://gcc.fyxm.net/releases/gcc-4.8.1/gcc-4.8.1.tar.bz2
# https://waf.googlecode.com/files/waf-1.7.11
# http://www.nic.funet.fi/pub/gnu/ftp.gnu.org/pub/gnu/gsl/gsl-1.15.tar.gz
# ftp://xmlsoft.org/libxml2/libxml2-2.9.1.tar.gz
# http://geuz.org/gmsh/src/gmsh-2.7.1-source.tgz
# http://www.cmake.org/files/v2.8/cmake-2.8.11.2.tar.gz
# http://apache-mirror.rbc.ru/pub/apache//apr/apr-1.4.8.tar.gz
# http://apache-mirror.rbc.ru/pub/apache//apr/apr-util-1.5.2.tar.gz
# http://www.sai.msu.su/apache/logging/log4cxx/0.10.0/apache-log4cxx-0.10.0.tar.gz
# http://www.vtk.org/files/release/5.2/vtk-5.2.1.tar.gz
# ftp://ftp.freedesktop.org/pub/mesa/older-versions/7.x/7.7.1/MesaLib-7.7.1.tar.gz
# http://python.org/ftp/python/2.7.5/Python-2.7.5.tar.bz2
# http://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v4.0&type=source&os=all&downloadFile=ParaView-v4.0.1-source.tgz


echo "Downloading dependencies"
if [[ -z "$GCM3D_DEPS_PATH" ]]; then
    # Downloading from Dropbox «mirror»
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/md5
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/mpfr-3.1.2.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/gmp-5.1.2.tar.bz2
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/mpc-1.0.1.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/cloog-0.18.0.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/isl-0.11.1.tar.bz2
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/gcc-4.8.1.tar.bz2
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/waf-1.7.11
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/gsl-1.15.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/libxml2-2.9.1.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/gmsh-2.7.1-source.tgz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/cmake-2.8.11.2.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/apr-1.4.8.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/apr-util-1.5.2.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/apache-log4cxx-0.10.0.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/vtk-5.2.1.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/MesaLib-7.7.1.tar.gz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/Python-2.7.5.tar.bz2
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/ParaView-v4.0.1-source.tgz
    wget -nv https://dl.dropboxusercontent.com/u/19548163/mipt/gcm3d-deps/patches.tar.gz
else
    cp -v "$GCM3D_DEPS_PATH"/* .
fi;

echo "Downloading latest gcm3d source code"
git clone https://github.com/avasyukov/gcm-3d.git

echo "Creating build script"
get_build_script > build.sh
chmod +x build.sh

cd "$DIR"

echo "Creating bootstrap script"
get_bootstrap_script > gcm3d-bootstrap.sh
echo "Packing source code and dependencies"
tar -cvz -C /tmp/gcm-3d-bundle . >> gcm3d-bootstrap.sh
chmod +x gcm3d-bootstrap.sh
echo "Done"

exit

# ----- attachment 1: script to build gcm3d with all dependencies
#!/bin/bash
# dependencies: gcc, g++, make, patch, gfortran (and a bunch of paraview dependencies)

die() {
    echo "$1" 2>&1
    exit 1
}

prepare_env() {
    echo "Preparing environment"
    if [[ -z "$GCM3D_INSTALL_PATH" ]]; then
        export GCM3D_INSTALL_PATH="$HOME/gcm3d"
    fi;
    export PATH="$GCM3D_INSTALL_PATH/bin:$PATH"
    export LD_LIBRARY_PATH="$GCM3D_INSTALL_PATH/lib:$LD_LIBRARY_PATH"
    export LD_LIBRARY_PATH="$GCM3D_INSTALL_PATH/lib64:$LD_LIBRARY_PATH"
    export LD_LIBRARY_PATH="$GCM3D_INSTALL_PATH/lib/vtk-5.2:$LD_LIBRARY_PATH"
}

unpack() {
    echo "Checking md5 sums"
    echo "Unpacking"
    {
        md5sum -c md5  || return -1           && \
        cd "$DIR"                             && \
        mkdir patches                         && \
        tar -xf patches.tar.gz -C patches     && \
        tar -xf mpfr-3.1.2.tar.gz             && \
        tar -xf gmp-5.1.2.tar.bz2             && \
        tar -xf mpc-1.0.1.tar.gz              && \
        tar -xf cloog-0.18.0.tar.gz           && \
        tar -xf isl-0.11.1.tar.bz2            && \
        tar -xf gcc-4.8.1.tar.bz2             && \
        tar -xf gsl-1.15.tar.gz               && \
        tar -xf libxml2-2.9.1.tar.gz          && \
        tar -xf cmake-2.8.11.2.tar.gz         && \
        tar -xf gmsh-2.7.1-source.tgz         && \
        tar -xf apr-1.4.8.tar.gz              && \
        tar -xf apr-util-1.5.2.tar.gz         && \
        tar -xf apache-log4cxx-0.10.0.tar.gz  && \
        tar -xf vtk-5.2.1.tar.gz              && \
        tar -xf MesaLib-7.7.1.tar.gz          && \
        tar -xf Python-2.7.5.tar.bz2          && \
        tar -xf ParaView-v4.0.1-source.tgz
    } &> "$DIR/unpacking.log"
}

install_gcc() {
    echo "Installing gcc"
    {
        mv "$DIR/gmp-5.1.2" "$DIR/gcc-4.8.1/gmp"      && \
        mv "$DIR/mpc-1.0.1" "$DIR/gcc-4.8.1/mpc"      && \
        mv "$DIR/mpfr-3.1.2" "$DIR/gcc-4.8.1/mpfr"    && \
        mv "$DIR/cloog-0.18.0" "$DIR/gcc-4.8.1/cloog" && \
        mv "$DIR/isl-0.11.1" "$DIR/gcc-4.8.1/isl"     && \
        chmod +x "$DIR/gcc-4.8.1/configure"           && \
        mkdir "$DIR/gccbuild"                         && \
        cd "$DIR/gccbuild"                            && \
        CXXFLAGS=-DHAVE_ICONV                            \
        ../gcc-4.8.1/configure                           \
            --prefix="$GCM3D_INSTALL_PATH"               \
            --disable-multilib                           \
            --enable-languages=c,c++                     \
            --enable-shared                              \
            --disable-werror                             \
            --disable-bootstrap                       && \
        make                                          && \
        make install
    } &> "$DIR/gcc.log"
}

install_waf() {
    echo "Installing waf"
    {
        cd "$DIR"                                           && \
        mkdir -p "$GCM3D_INSTALL_PATH/bin"                  && \
        patch -p0 -i "$DIR/patches/waf_python_binary.patch" && \
        mv waf-1.7.11 "$GCM3D_INSTALL_PATH/bin/waf"         && \
        chmod +x "$GCM3D_INSTALL_PATH/bin/waf"
    } &> "$DIR/waf.log"
}

install_gsl() {
    echo "Installing gsl"
    {
        cd "$DIR/gsl-1.15"                 && \
        ./configure                           \
            --prefix="$GCM3D_INSTALL_PATH" && \
        make                               && \
        make install
    } &> "$DIR/gsl.log"
}

install_libxml2() {
    echo "Installing libxml2"
    {
        cd "$DIR/libxml2-2.9.1"           && \
        ./configure                          \
            --prefix="$GCM3D_INSTALL_PATH"   \
            --without-python              && \
        make                              && \
        make install
    } &> "$DIR/libxml2.log"
}

install_cmake(){
    echo "Installing cmake"
    {
        cd "$DIR/cmake-2.8.11.2"           && \
        ./configure                           \
            --prefix="$GCM3D_INSTALL_PATH" && \
        make                               && \
        make install
    } &> "$DIR/cmake.log"
}

install_gmsh() {
    echo "Installing gmsh"
    {
        cd "$DIR/gmsh-2.7.1-source"                           && \
        patch                                                    \
            -p1                                                  \
            -i "$DIR/patches/gmsh-addruledfaces-return.patch" && \
        cmake                                                    \
            -DCMAKE_INSTALL_PREFIX="$GCM3D_INSTALL_PATH"         \
            -DENABLE_BUILD_DYNAMIC=1                             \
            -DENABLE_BUILD_SHARED=1                           && \
        make                                                  && \
        make install
    } &> "$DIR/gmsh.log"
}

install_apr() {
    echo "Installing apr"
    {
        cd "$DIR/apr-1.4.8"                && \
        ./configure                           \
            --prefix="$GCM3D_INSTALL_PATH" && \
        make                               && \
        make install
    } &> "$DIR/apr.log"
}

install_apr_util() {
    echo "Installing apr-util"
    {
        cd "$DIR/apr-util-1.5.2"                              && \
        ./configure                                              \
            --prefix="$GCM3D_INSTALL_PATH"                       \
            --with-apr="$GCM3D_INSTALL_PATH/bin/apr-1-config" && \
        make                                                  && \
        make install
    } &> "$DIR/apr_util.log"
}

install_log4cxx() {
    echo "Installing log4cxx"
    {
        cd "$DIR/apache-log4cxx-0.10.0"                            && \
        patch                                                         \
            -p1                                                       \
            -i "$DIR/patches/cppFolder_stringInclude.patch"        && \
        patch                                                         \
            -p1                                                       \
            -i "$DIR/patches/exampleFolder_stringInclude.patch"    && \
        ./configure                                                   \
            --prefix="$GCM3D_INSTALL_PATH"                            \
            --with-apr="$GCM3D_INSTALL_PATH/bin/apr-1-config"         \
            --with-apr-util="$GCM3D_INSTALL_PATH/bin/apu-1-config" && \
        make                                                       && \
        make install
    } &> "$DIR/log4cxx.log"
}

install_vtk() {
    echo "Installing vtk"
    {
        cd "$DIR/VTK"                                    && \
        cmake                                               \
            -DCMAKE_INSTALL_PREFIX="$GCM3D_INSTALL_PATH"    \
            -DVTK_USE_RENDERING=0                           \
            -DVTK_USE_TK=0                                  \
            -DVTK_WRAP_PYTHON=0                             \
            -DVTK_WRAP_TCL=0                                \
            -DVTK_WRAP_JAVA=0                               \
            -DBUILD_SHARED_LIBS=1                           \
            -DVTK_INSTALL_NO_DEVELOPMENT=0               && \
        make clean                                       && \
        make                                             && \
        make install
    } &> "$DIR/vtk.log"
}

install_gcm3d() {
    echo "Installing gcm3d"
    {
        cd "$DIR/gcm-3d"                                        && \
        CXX="$GCM3D_INSTALL_PATH/bin/g++"                          \
        waf configure                                              \
            --prefix="$GCM3D_INSTALL_PATH"                         \
            --without-tests                                        \
            --without-headers                                      \
            --without-resources                                    \
            --includepath="$GCM3D_INSTALL_PATH/include/libxml2"    \
            --includepath="$GCM3D_INSTALL_PATH/include/vtk-5.2"    \
            --includepath="$GCM3D_INSTALL_PATH/include/gmsh"       \
            --includepath="$GCM3D_INSTALL_PATH/include"            \
            --libpath="$GCM3D_INSTALL_PATH/lib/vtk-5.2"            \
            --libpath="$GCM3D_INSTALL_PATH/lib"                    \
            --libpath="$GCM3D_INSTALL_PATH/lib64"               && \
        waf build                                               && \
        waf install
    } &> "$DIR/gcm3d.log"
}

install_mesa() {
    echo "Installing mesa3D"
    {
        cd "$DIR/Mesa-7.7.1"               && \
        ./configure                           \
            --prefix="$GCM3D_INSTALL_PATH"    \
            --with-driver=xlib                \
            --enable-gl-osmesa                \
            --disable-glw                  && \
        make                               && \
        make install
    } &> "$DIR/mesa3d.log"
}

install_python() {
    echo "Installing python"
    {
        cd "$DIR/Python-2.7.5"             && \
        ./configure                           \
            --prefix="$GCM3D_INSTALL_PATH"    \
            --enable-shared                && \
        make                               && \
        make install

    } &> "$DIR/python.log"
}

install_paraview() {
    echo "Installing paraview"
    {
        mkdir "$DIR/pvbuild"                                       && \
        cd "$DIR/pvbuild"                                          && \
        CXX=mpicxx                                                    \
        CC=mpicc                                                      \
        cmake                                                         \
            -DCMAKE_INSTALL_PREFIX="$GCM3D_INSTALL_PATH"              \
            -DPARAVIEW_BUILD_QT_GUI=0                                 \
            -DPARAVIEW_USE_MPI=1                                      \
            -DBUILD_SHARED_LIBS=1                                     \
            -DPARAVIEW_ENABLE_PYTHON=1                                \
            -DVTK_USE_X=0                                             \
            -DCMAKE_BUILD_TYPE=Release                                \
            -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=0                    \
            -DBUILD_TESTING=0                                         \
            -DVTK_OPENGL_HAS_OSMESA=1                                 \
            -DOSMESA_INCLUDE_DIR="$GCM3D_INSTALL_PATH/include"        \
            -DOSMESA_LIBRARY="$GCM3D_INSTALL_PATH/lib/libOSMesa.so"   \
            -DOPENGL_INCLUDE_DIR="$GCM3D_INSTALL_PATH/include"        \
            -DOPENGL_gl_LIBRARY="$GCM3D_INSTALL_PATH/lib/libGL.so"    \
            -DOPENGL_glu_LIBRARY="$GCM3D_INSTALL_PATH/lib/libGLU.so"  \
            "$DIR/ParaView-v4.0.1-source"                          && \
        make                                                       && \
        make install
    } &> "$DIR/paraview.log"
}

patch_bashrc() {
    if [[ ! -f "$HOME/.bashrc" ]] || ! grep -q gcm3d_bashrc "$HOME/.bashrc"; then
        echo 'source "$HOME"/.gcm3d_bashrc' > "$HOME/.bashrc"
    fi
    echo "export GCM3D_INSTALL_PATH=\"$GCM3D_INSTALL_PATH\""                         >  "$HOME/.gcm3d_bashrc"
    echo 'export PATH="$GCM3D_INSTALL_PATH/bin:$PATH"'                               >> "$HOME/.gcm3d_bashrc"
    echo 'export LD_LIBRARY_PATH="$GCM3D_INSTALL_PATH/lib:$LD_LIBRARY_PATH"'         >> "$HOME/.gcm3d_bashrc"
    echo 'export LD_LIBRARY_PATH="$GCM3D_INSTALL_PATH/lib64:$LD_LIBRARY_PATH"'       >> "$HOME/.gcm3d_bashrc"
    echo 'export LD_LIBRARY_PATH="$GCM3D_INSTALL_PATH/lib/vtk-5.2:$LD_LIBRARY_PATH"' >> "$HOME/.gcm3d_bashrc"
}

install() {
    SCRIPT="`readlink -f $0`" \
    DIR="`dirname "$SCRIPT"`" \
    prepare_env            && \
    unpack                 && \
    install_gcc            && \
    install_waf            && \
    install_cmake          && \
    install_apr            && \
    install_apr_util       && \
    install_log4cxx        && \
    install_gsl            && \
    install_libxml2        && \
    install_vtk            && \
    install_gmsh           && \
    install_gcm3d          && \
    install_mesa           && \
    install_python         && \
    install_paraview       && \
    patch_bashrc           || die "Build failed. See corresponding log file for details."
}

clean() {
    echo "Removing unneeded source code"
    {
    rm -rf                            \
        "$DIR/gcc-4.8.1"              \
        "$DIR/gccbuild"               \
        "$DIR/gsl-1.15"               \
        "$DIR/libxml2-2.9.1"          \
        "$DIR/cmake-2.8.11.2"         \
        "$DIR/gmsh-2.7.1-source"      \
        "$DIR/apr-1.4.8"              \
        "$DIR/apr-util-1.5.2"         \
        "$DIR/apache-log4cxx-0.10.0"  \
        "$DIR/VTK"                    \
        "$DIR/Mesa-7.7.1"             \
        "$DIR/Python-2.7.5"           \
        "$DIR/ParaView-v4.0.1-source" \
        "$DIR/pvbuild"
    } &> clean.log
}

if [[ "x$1" = "xinstall" ]]; then
    install
fi

exit

# ----- attachment 2: bootstrap script
#!/bin/bash
SCRIPT=`readlink -f $0`

if [[ -z "$GCM3D_INSTALL_TMP" ]]; then
    export GCM3D_INSTALL_TMP="$HOME/tmp"
fi

echo "Preparing installation..."

mkdir -p "$GCM3D_INSTALL_TMP"
cd "$GCM3D_INSTALL_TMP"

rm -rf gcm-3d-bootstrap && mkdir gcm-3d-bootstrap && cd gcm-3d-bootstrap
sed -e '1,/^exit$/d' "$SCRIPT" | tar xzf - && ./build.sh install

exit
