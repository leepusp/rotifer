-- Rotifer's shell configuration

-- ------------------------------
-- Lua Funcions
-- ------------------------------
local function read_file(path)
    local file = open(path, "rb") -- r read mode and b binary mode
    if not file then return nil end
    local content = file:read "*a" -- *a or *all reads the whole file
    file:close()
    return content
end

-- ------------------------------
-- Local variables
-- ------------------------------
local myversion = "0.3"
local myfn  = myFileName()
local mymfn = myModuleFullName()
local mydir = "none"
if realpath ~= nil then
	mydir = realpath(dirname(realpath(myfn)) .. "/../../../../")
else
	mydir = abspath(dirname(abspath(myfn)) .. "/../../../../")
end
local datadir = mydir .. "/share/rotifer/data"

-- ------------------------------
-- Environment variables
-- ------------------------------
prepend_path("PYTHONPATH", mydir .. "/lib")
prepend_path("PERL5LIB", mydir .. "/perl/lib")
if (os.getenv("DATABASES") ~= nil) then
	datadir = os.getenv("DATABASES")
end
setenv("ROTIFER_DATA",datadir)

-- PATH
prepend_path("PATH", mydir .. "/bin")

-- ------------------------------
-- Shell Functions
-- ------------------------------

-- acc2pfam
set_shell_function("acc2pfam"," \
    SEQ=$(efetch -db protein -format fasta -id $1);\
    ( echo \"$SEQ\" | hmmscan --cpu 4 $ROTIFER_DATA/pfam/Pfam-A.hmm - | hmmer2table -c model=version;\
    echo \"$SEQ\" | phobius 2> /dev/null | phobius2table -e 0.0101 ) | domain2architecture -e 0.0101 | architecture2table | padtable\
","")

-- acc2profiledb
set_shell_function("acc2profiledb"," \
    ( efetch -db protein -format fasta -id $1 | rpsblast -db $ROTIFER_DATA/rpsdb/allprofiles | blast2table -c profiledb -s;\
    efetch -db protein -format fasta -id $1 | phobius 2> /dev/null | phobius2table -e 0.0101 ) | domain2architecture -e 0.0101 | architecture2table | padtable\
","")

-- aln_order_by_tree
set_shell_function("aln_order_by_tree"," \
    aln=$1;\
    tree=$2;\
    tmp=$(mktemp -p /tmp aln_order_by_tree.XXXXXX);\
    treeutil -of leaves $tree > $tmp;\
    aln2seqrows -r '\\t' $aln | tjoin -r '\\n' -i1 1 -i2 1 -1 1 -f1 '$F[1] ne \"ID\"' -f '$F[0]=\">$F[0]\";1' -p $tmp -;\
    rm -f $tmp\
","")

-- distinct_dom_per_prot
set_shell_function("distinct_dom_per_prot"," \
    arch=$1;\
    shift;\
    architecture2table $arch | tgroup -g 0 -a 1=count:distinct -i 0..1 \"$@\"\
","")

-- dom_count
set_shell_function("dom_count"," \
    arch=$1;\
    shift;\
    architecture2table $arch | tgroup -f '$F[2]=$F[1];1' -g 0 -g 1 -a 2=count -i 0..2 \"$@\"\
","")

-- dom_per_prot
set_shell_function("dom_per_prot"," \
    arch=$1;\
    shift;\
    architecture2table $arch | tgroup -g 0 -a 1=count -i 0..1 \"$@\"\
","")

-- tdesc
set_shell_function("tdesc"," \
    sep=$2;\
    if [ \"$sep\" == \"\" ]; then\
        sep='\\t';\
    fi;\
    head -n1 $1 | ttranspose -s $sep | tfilter -f 'unshift(@F,$H{c}++);1'\
","")

-- ------------------------------
-- Whatis
-- ------------------------------
whatis("Module: " .. mymfn)
whatis("Version: " .. myversion)
whatis("Filename: " .. myfn)
whatis("Directory: " .. mydir)
whatis("Data directory: " .. datadir)
whatis("Category: data analysis ")
whatis("Description: Rotifer's shell setup.")
whatis("URL: https://github.com/leepusp/rotifer ")

-- Help message
local helpMsg = [[
ROTIFER
=======

Rapid Open-source Tools and Infrastructure For data Exploration and Research
----------------------------------------------------------------------------

ROTIFER is a multi-language collection of high-level programming libraries
for the development of data analysis pipelines, mostly targeting problems in
comparative genomics and the computational analysis of biological sequences.

It also provides a collection of easy to use command line tools based on this framework.

Current installation directory: ]] .. mydir

