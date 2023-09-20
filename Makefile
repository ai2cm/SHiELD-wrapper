VERSION ?= $(shell git rev-parse HEAD)

build:
	docker build -f docker/Dockerfile -t shield-wrapper:$(VERSION) .

test:
	docker/docker_run.sh \
		-v $(shell pwd)/wrapper/tests:/wrapper/tests \
		shield-wrapper:$(VERSION) make -C /wrapper test

test_regtest_reset:
	docker/docker_run.sh \
		-v $(shell pwd)/wrapper/tests:/wrapper/tests \
		shield-wrapper:$(VERSION) make -C /wrapper test_regtest_reset

enter:
	docker/docker_run.sh \
		-v $(shell pwd)/wrapper/tests:/wrapper/tests \
		-it shield-wrapper:$(VERSION) bash

lock_pip:
	pip-compile requirements.in
