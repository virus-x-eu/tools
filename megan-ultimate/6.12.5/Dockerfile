FROM ubuntu:rolling
ENV MEGAN_VERSION=6_12_5 PATH="/usr/local/megan/tools:${PATH}" DEBIAN_FRONTEND=noninteractive
ADD MEGAN-headless /usr/local/bin/MEGAN-headless
RUN apt-get update && \
    apt-get install -yq xvfb libgtk2.0-0 openjfx wget && \
    wget --progress=dot:mega "http://ab.inf.uni-tuebingen.de/data/software/megan6/download/MEGAN_Ultimate_unix_${MEGAN_VERSION}.sh" -O /tmp/megan_installer.sh && \
    echo 'o\n1\n/usr/local/megan\n1,3\ny\n/usr/local/bin\n5\n12000\nn' | sh /tmp/megan_installer.sh && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*