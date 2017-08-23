FROM java:8-jdk-alpine
RUN apk add --update wget bash && \
    cd /usr/local/share/ && \
    wget "https://bibiserv.cebitec.uni-bielefeld.de/resources/emgb2/EMGBDataWarehouse-0.3.0.jar" && \
    apk del wget bash && \
    rm -rf /var/cache/apk/*
