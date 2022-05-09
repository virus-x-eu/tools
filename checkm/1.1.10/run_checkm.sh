#!/bin/bash

echo "{\"dataRoot\": \"/srv/whitlam/bio/db/checkm_data/1.0.0\", \"remoteManifestURL\": \"https://data.ace.uq.edu.au/public/CheckM_databases/\", \"manifestType\": \"CheckM\", \"remoteManifestName\": \".dmanifest\", \"localManifestName\": \".dmanifest\"}" > /tmp/DATA_CONFIG
checkm data setRoot "${HOME}/checkm_data" "/tmp/DATA_CONFIG"
checkm ${@}