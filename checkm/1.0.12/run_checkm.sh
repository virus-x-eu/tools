#!/bin/bash

echo "${HOME}/checkm_data_v${CHECKM_VERSION}/" | checkm data setRoot "${HOME}/checkm_data_v${CHECKM_VERSION}/"
checkm ${@}