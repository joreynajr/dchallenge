# setting important paths
hichip_db_py="/mnt/BioAdHoc/Groups/vd-ay/jreyna/software/mambaforge/envs/hichip-db/bin/python"
hichip_db_R="/mnt/BioApps/R/3.6.1/bin/R"
hichip_db_Rscript="/mnt/BioApps/R/3.6.1/bin/Rscript"
bgzip="/mnt/BioApps/tabix/tabix-0.2.6/bgzip"
tabix="/mnt/BioApps/tabix/tabix-0.2.6/tabix"

#Setting helper functions
function path_exists(){
    if [[ -e "$1" ]]
    then
        echo "Is-Present"
    else
        echo "Not-Present"
    fi
}
