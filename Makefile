VERSION ?= $(shell git rev-parse HEAD)

build:
	docker build -f docker/Dockerfile -t shield-wrapper:$(VERSION) .
