#! /bin/bash -x

type=${1:-clean}

if [ ${type} = "clean" ]; then
   cd 1x96
   rm -rf ./scale3_*.cnf ./k_stgpjsub_*.sh
   cd -

   cd 16x96
   rm -rf ./scale3_*.cnf ./k_stgpjsub_*.sh
   cd -

   cd 32x96
   rm -rf ./scale3_*.cnf ./k_stgpjsub_*.sh
   cd -

   cd 64x96
   rm -rf ./scale3_*.cnf ./k_stgpjsub_*.sh
   cd -

   cd 128x96
   rm -rf ./scale3_*.cnf ./k_stgpjsub_*.sh
   cd -

   cd 192x96
   rm -rf ./scale3_*.cnf ./k_stgpjsub_*.sh
   cd -

   cd 256x96
   rm -rf ./scale3_*.cnf ./k_stgpjsub_*.sh
   cd -
elif [ ${type} = "new" ]; then
   cd 1x96
   sh convert_nodenum.sh
   cd -

   cd 16x96
   sh convert_nodenum.sh
   cd -

   cd 32x96
   sh convert_nodenum.sh
   cd -

   cd 64x96
   sh convert_nodenum.sh
   cd -

   cd 128x96
   sh convert_nodenum.sh
   cd -

   cd 192x96
   sh convert_nodenum.sh
   cd -

   cd 256x96
   sh convert_nodenum.sh
   cd -
fi