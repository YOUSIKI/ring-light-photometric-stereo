FROM ghcr.io/gnu-octave/octave:latest

RUN octave -q --eval "pkg install -forge image"
RUN octave -q --eval "pkg install -forge control"
RUN octave -q --eval "pkg install -forge signal"
