# setting important paths
hichip_db_py="/mnt/BioAdHoc/Groups/vd-ay/jreyna/software/mambaforge/envs/hichip-db/bin/python"

#Setting helper functions
function path_exists(){
    if [[ -e "$1" ]]
    then
        echo "Is-Present"
    else
        echo "Not-Present"
    fi
}
