mkdir -p results/FARM/inputs/x509_user_proxy; voms-proxy-init --voms cms:/cms/becms  -valid 192:00 --out results/FARM/inputs/x509_proxy
export X509_USER_PROXY=results/FARM/inputs/x509_proxy;voms-proxy-init --noregen;
