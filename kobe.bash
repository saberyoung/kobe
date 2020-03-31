#!/bin/bash

clear
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo "... KOBE package initialize @${DIR} ..."

export PYTHONPATH="${DIR}/src/"
