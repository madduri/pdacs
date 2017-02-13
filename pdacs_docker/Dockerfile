FROM quay.io/bgruening/galaxy:gcc2016

MAINTAINER madduri@anl.gov

ENV GALAXY_CONFIG_BRAND=PDACS \
    GALAXY_LOGGING=full

RUN apt-get -y install sqlite3 libsqlite3-dev

WORKDIR /galaxy-central

ADD tool_conf.xml /galaxy-central/tool_conf.xml

RUN rm -rf tools

COPY tools /galaxy-central/tools

ENV GALAXY_CONFIG_TOOL_CONFIG_FILE /galaxy-central/tool_conf.xml
# overwrite current welcome page
ADD welcome.html $GALAXY_CONFIG_DIR/web/welcome.html
ADD Slide1.jpg $GALAXY_CONFIG_DIR/web/Slide1.jpg

# Mark folders as imported from the host.
VOLUME ["/export/", "/var/lib/docker"]

