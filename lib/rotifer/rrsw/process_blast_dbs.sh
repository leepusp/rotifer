#!/bin/bash

### START OF TEMPLATE HEADER: do not change this section! 
# Input for run_or_show must avoid interpolation of pipes and redirections!
function run_or_show() {
	if [ "$IS_TEST" == "1" ]; then echo "$@"; else eval "$@"; fi
}
export IS_TEST=$1 SOURCE_DIR=$2 DEST_DIR=$3 IWD=`pwd`
if [ "$IS_TEST" == "1" ]; then cd $DEST_DIR; fi
if [ -f /home/linuxbrew/.linuxbrew.sh ]; then source /home/linuxbrew/.linuxbrew.sh; fi
### END OF TEMPLATE HEADER

### START OF CUSTOMIZABLE SECTION: you may change or add rows below
base=$(dirname $DEST_DIR)
if [ "$IS_TEST" == "1" ]; then echo command $0, parameters: $IS_TEST $SOURCE_DIR $DEST_DIR, PWD: $(pwd), base: $base; fi
if [ ! -d "$base/fadb"   ]; then run_or_show mkdir -p "$base/fadb"; fi
for f in $(\ls -1 db/FASTA/* 2> /dev/null | grep -Fv 'db/FASTA/*' | grep -v '\.md5$' | grep -v '\.ssi$' | grep -v '\.p..$')
do
	target=$(basename $f)
	if [ "$target" == "nt" ]; then DBTYPE="nucl"; fi
	run_or_show "sed -i '/^>/!s/-/X/g' \"$f\""
	run_or_show esl-sfetch --index $f
done
if [   -d "$base/freeze" ]; then run_or_show  ln -s $base/freeze .; fi
### END OF CUSTOMIZABLE SECTION: your code ends here!

### START OF TEMPLATE FOOT: do not change this section! 
cd $IWD
exit $?
### END OF TEMPLATE FOOT: any rows below will have no effect!
