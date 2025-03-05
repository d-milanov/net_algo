#!/bin/bash

echo Processing configured in the conf/sample_conf.txt file takes 1-15 minutes depending on your computer.
echo You can track progress in the BF_dbg and Net_dbg files.
echo Program execution results are in the results/Rho_02_07_1 file.

java  -Xms2548M -Xmx2548M -cp app/build/libs/app.jar  celestial.orbits.net_algo.Test conf/sample_conf.txt  >results/Rho_02_07_1
