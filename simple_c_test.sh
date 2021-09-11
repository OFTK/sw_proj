red='\033[0;31m'
blue='\033[0;34m'

echo "${blue}\nWAM TEST RESULTS:"
./spkmeans 0 wam tests/input.csv > tests/our_wam_output.csv
echo "${red}"
diff tests/rami_wam_output.csv tests/our_wam_output.csv

echo "${blue}\nDDG TEST RESULTS:"
./spkmeans 0 ddg tests/input.csv > tests/our_ddg_output.csv
echo "${red}"
diff tests/rami_ddg_output.csv tests/our_ddg_output.csv

echo "${blue}\nLNORM TEST RESULTS:"
./spkmeans 0 lnorm tests/input.csv > tests/our_lnorm_output.csv
echo "${red}"
diff tests/rami_lnorm_output.csv tests/our_lnorm_output.csv

echo "${blue}\nSPK TEST RESULTS:"
./spkmeans 0 spk tests/input.csv > tests/our_spk_c_output.csv
echo "${red}"
diff tests/rami_spk_c_output.csv tests/our_spk_c_output.csv

echo "${blue}\nJACOBI TEST RESULTS:"
./spkmeans 0 jacobi tests/jacobi_input.csv > tests/our_jacobi_output.csv
echo "${red}"
diff tests/rami_jacobi_output.csv tests/our_jacobi_output.csv
