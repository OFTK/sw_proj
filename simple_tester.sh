red='\033[0;31m'
blue='\033[0;34m'

rm tests/our_wam_output.csv
rm tests/our_ddg_output.csv
rm tests/our_lnorm_output.csv
rm tests/our_jacobi_output.csv
rm tests/our_spk_c_output.csv

echo "${blue}\nWAM TEST RESULTS:"
./spkmeans 0 wam tests/accurate_rami_inputs > tests/our_wam_output.csv
echo "${red}"
diff tests/rami_wam_output.csv tests/our_wam_output.csv

echo "${blue}\nDDG TEST RESULTS:"
./spkmeans 0 ddg tests/accurate_rami_inputs > tests/our_ddg_output.csv
echo "${red}"
diff tests/rami_ddg_output.csv tests/our_ddg_output.csv

echo "${blue}\nLNORM TEST RESULTS:"
./spkmeans 0 lnorm tests/accurate_rami_inputs > tests/our_lnorm_output.csv
echo "${red}"
diff tests/rami_lnorm_output.csv tests/our_lnorm_output.csv

echo "${blue}\nSPK TEST RESULTS:"
./spkmeans 0 spk tests/accurate_rami_inputs > tests/our_spk_c_output.csv
echo "${red}"
diff tests/rami_spk_c_output.csv tests/our_spk_c_output.csv

echo "${blue}\nJACOBI TEST RESULTS:"
./spkmeans 0 jacobi tests/jacobi_input.csv > tests/our_jacobi_output.csv
echo "${red}"
diff tests/rami_jacobi_output.csv tests/our_jacobi_output.csv
