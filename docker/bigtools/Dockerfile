FROM rust:1.93.0

WORKDIR /opt
RUN git clone https://github.com/jackh726/bigtools.git \
	&& cd bigtools \
	&& git checkout "3844b58" \
	&& cargo build --release

	ENV PATH "${PATH}:/opt/bigtools/target/release"