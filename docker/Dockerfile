FROM centos:6.10

COPY CentOS-Base.repo /etc/yum.repos.d/


RUN yum update -y && \
    yum upgrade -y && \
    yum install -y xz git libzip-devel curl wget zlib-devel xz-devel bzip2-devel openssl-devel libcurl-devel

RUN wget --no-check-certificate http://people.centos.org/tru/devtools-2/devtools-2.repo -O /etc/yum.repos.d/devtools-2.repo && \
    yum install -y devtoolset-2-gcc

RUN yum install -y devtoolset-2-binutils devtoolset-2-gcc-c++ && \
    source scl_source enable devtoolset-2 && \
    echo "source scl_source enable devtoolset-2" >> ~/.bashrc && \
    echo "source scl_source enable devtoolset-2" >> ~/.bash_profile


RUN source scl_source enable devtoolset-2 && \
    wget --quiet --no-check-certificate https://ftp.gnu.org/gnu/m4/m4-1.4.18.tar.gz && \
    tar xzf m4-1.4.18.tar.gz && cd m4* && ./configure && make && make install && cd .. && rm -rf m4* && \
    wget --no-check-certificate --quiet http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz && \
    tar xzf autoconf-2.69.tar.gz && \
    cd autoconf* && ./configure && make && make install && cd .. && rm -rf autoconf* && \
    git clone --depth 1 https://github.com/ebiggers/libdeflate.git && \
    cd libdeflate && make -j 2 CFLAGS='-fPIC -O3' libdeflate.a && \
    cp libdeflate.a /usr/local/lib && cp libdeflate.h /usr/local/include && \
    cd .. && rm -rf libdeflate && \
    git clone --recursive https://github.com/samtools/htslib && \
    cd htslib && git checkout 1.13 && autoheader && autoconf && \
    ./configure --disable-libcurl --with-libdeflate && \
    cd .. && make -j4 CFLAGS="-fPIC -O3" -C htslib install && \
    git clone https://github.com/samtools/bcftools && \
    cd bcftools && git checkout 1.13 && autoheader && autoconf && \
    ./configure --disable-libcurl --with-libdeflate --disable-bcftools-plugins && \
    cd .. && make -j4 PLUGINS_ENABLED=no CFLAGS="-O3" -C bcftools install && \
    echo "/usr/local/lib" >> etc/ld.so.conf


RUN cd / \
    && wget --no-check-certificate --progress=dot:giga https://github.com/nim-lang/nightlies/releases/download/2020-12-22-devel-297c8e403d110dd872e070563328f4e0c734cd01/nim-1.5.1-linux_x64.tar.xz \
    && tar Jxf nim-1.5.1-linux_x64.tar.xz \
    && rm -f nim-1.5.1-linux_x64.tar.xz \
    && mv nim-1.5.1 nim && \
    echo 'PATH=/nim/bin:/root/.nimble/bin:$PATH' >> ~/.bashrc && \
    echo 'PATH=/nim/bin:/root/.nimble/bin:$PATH' >> ~/.bash_profile


RUN cd / && \
    ldconfig && \
    source scl_source enable devtoolset-2 && \
    export PATH=/nim/bin:$PATH && \
    git clone https://github.com/iffy/nim-argparse && cd nim-argparse && git checkout v0.10.1 && nimble install -y && \
    cd .. && rm -rf nim-argparse && \
    git clone --depth 1 git://github.com/brentp/slivar.git && \
    cd slivar && \
    sed -i 's/-no-pie//' nim.cfg && \
    nimble install -y &&  \
    nim c -d:release -d:danger -o:/usr/bin/slivar src/slivar && \
    mkdir -p /opt/slivar/ && cp js/slivar-functions.js /opt/slivar/ && \
    rm -rf /slivar /duktape-nim && slivar expr -h
