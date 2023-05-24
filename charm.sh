while getopts i:S: flag
do
    case "${flag}" in
        i) ini=${OPTARG};;
        S) stage=${OPTARG};;
    esac
done

python3 scripts/charm_manager.py -S $stage -i $ini