#!/bin/bash

# Taken largely from fv3net:
# https://github.com/ai2cm/fv3net/blob/master/tools/docker-run

set -e
dockerArgs=()

# Google authentication
if [[ -n $GOOGLE_APPLICATION_CREDENTIALS ]] && \
   [[ "$(jq -r .type $GOOGLE_APPLICATION_CREDENTIALS)" == "service_account" ]]
then
    >&2 echo "Using service account credentials"
    dockerArgs+=(-v "${GOOGLE_APPLICATION_CREDENTIALS}:/tmp/key.json")
    dockerArgs+=(-e "GOOGLE_APPLICATION_CREDENTIALS=/tmp/key.json")
else
    >&2 echo "Using user credentials"
    dockerArgs=(-v ~/.config/gcloud:/root/.config/gcloud)
fi

docker run -e FSSPEC_GS_REQUESTER_PAYS=vcm-ml "${dockerArgs[@]}" $@
