build-linux:
	docker run --rm --user ${uid -u}:${uid -g} \
         -v ${PWD}:/usr/src/preqc-pack \
         -w /usr/src/preqc-pack \
         musl-build cargo build --release

build-docker:
	docker build -t musl-build .