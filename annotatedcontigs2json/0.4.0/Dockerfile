FROM golang:latest as builder
WORKDIR /tmp/
COPY annotatedcontigs2json.go /tmp/
RUN go get github.com/biogo/hts/bam && \
    go get github.com/shenwei356/bio/seqio/fastx && \
    GOARCH=amd64 GOOS=linux go build -ldflags "-linkmode external -extldflags -static -w" /tmp/annotatedcontigs2json.go

FROM alpine
COPY --from=builder /tmp/annotatedcontigs2json /usr/local/bin/