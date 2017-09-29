#!/bin/bash
echo "${1}" | checkm data setRoot "${1}"
checkm ${@:2}