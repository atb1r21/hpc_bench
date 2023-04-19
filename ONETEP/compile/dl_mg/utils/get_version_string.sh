branchname=$(git status -z -b -uno | awk '{print $2}')
dl_mg_src=$(git rev-parse --show-toplevel)/src
# regexp that should catch only the version string
re='^ *char.*parameter.*version_string *='
# check if their is some anomaly
nlines=$(grep -c -i "$re" "$dl_mg_src"/dl_mg_info.F90)
if [ ! x"$nlines" = x1 ]
then
    echo "got 0 or more than 1 line with version_string in dl_mg_info.F90. Please investigate!"
    grep -i -n "$re" "$dl_mg_src"/dl_mg_info.F90
    exit 1
fi

vs_line="$(grep -i "$re" "$dl_mg_src"/dl_mg_info.F90)"
#echo "vs line $vs_line"

# strip the begining of the version line and a possible comment
curent_version=$(echo "$vs_line" | sed -e 's/.*version_string=\(.*\)/\1/' -e 's/"//g' -e 's/[ (].*//' -e 's/!.*//')

current_version_with_date=$(echo "$vs_line"| sed -e 's/.*version_string=\(.*\)/\1/' -e 's/"//g' -e 's/!.*//')



# test if the version string was updated
# useful for pre-commit
if [ "$1" = -t ]
then
    if [ "$branchname" = releases ]
    then
	head_version=$(git show "$branchname":src/dl_mg_info.F90 | grep -i "$re" | sed -e 's/.*version_string=\(.*\)/\1/' -e 's/"//g' -e 's/[ (].*//' -e 's/!.*//')
	if [ "$curent_version" = "$head_version" ]
	then
	    echo "Error: I guess that you forgot to update the version string after working on branch releases"
	    exit 2
# test script will be called below
	fi
    fi
else
    case $1 in
    -with-date) echo "$current_version_with_date" ;;
        # alternative date from git
        # $(git log -n1 --format="(%cd)")"
    -with-hash)
	# add git hash to version string
	echo "$curent_version"-$(git log -n1 --format="%h (%cd)") ;;
    *) echo "$curent_version" ;;
    esac
fi
