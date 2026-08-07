@echo off
echo Processing configured in the conf/sample_conf.yaml file takes 1-15 minutes depending on your computer.
echo Program execution results are in the results/Rho_02_07_1 file.
@echo on
java  -Xms2548M -Xmx2548M -cp build/libs/net_algo-1.0-SNAPSHOT.jar  net.exmachine.app.orbits.performance.Test conf/timing/02_Rho_07.yaml > results/Rho_02_07_1
rem java  -Xms2548M -Xmx2548M -cp build/libs/net_algo-1.0-SNAPSHOT.jar  net.exmachine.app.orbits.performance.Test conf/timing/01_D_H_04.yaml >results/01_D_H_04
rem java  -Xms2548M -Xmx2548M -cp build/libs/net_algo-1.0-SNAPSHOT.jar  net.exmachine.app.orbits.performance.Test conf/timing/01_D2_05.yaml >results/01_D2_05
rem java  -Xms2548M -Xmx2548M -cp build/libs/net_algo-1.0-SNAPSHOT.jar  net.exmachine.app.orbits.performance.Test conf/timing/01_Rho_07.yaml >results/01_Rho_07
rem java  -Xms2548M -Xmx2548M -cp build/libs/net_algo-1.0-SNAPSHOT.jar  net.exmachine.app.orbits.performance.Test conf/timing/conf_01_D2_06_mpcorb.yaml >results/conf_01_D2_06_mpcorb
