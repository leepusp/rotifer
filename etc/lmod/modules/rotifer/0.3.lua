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
if (os.getenv("DATABASES") ~= nil) then
	datadir = os.getenv("DATABASES")
end
setenv("ROTIFER_DATA",datadir)

-- PATH
-- acc2pfam, acc2profiledb, aln_order_by_tree, distinct_dom_per_prot,
-- dom_count, dom_per_prot and tdesc, formerly shell functions defined
-- here, are commands in bin/
prepend_path("PATH", mydir .. "/bin")

-- Migration notice (ROTIFER's Perl tools were replaced by Python)
if (mode() == "load") then
	LmodMessage("ROTIFER: the Perl command line tools were replaced by Python. Removed tools print what to use instead.")
end

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

