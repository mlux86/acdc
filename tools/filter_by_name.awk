# filter_by_name.awk
# Argument 1: Name of file containing list of target contig identifiers (one per line, without leading '>' symbol)
# Argument 2: Name of fasta file to filter
# Output: Contigs from fasta file that match any of target contig identifiers

# trim functions
function ltrim(s) { sub(/^[ \t\r\n]+/, "", s); return s }
function rtrim(s) { sub(/[ \t\r\n]+$/, "", s); return s }
function trim(s) { return rtrim(ltrim(s)); }

# populate target contig names (keys) from first file
NR == FNR {
    keys[">"trim($0)]
    next
}

# for each contig name in second file
# set keep = 1 if it is present in keys variable
# and keep = 0 if not
$0 ~ /^>/ {
    contig = trim($0)
    if (contig in keys)
        keep = 1
    else
        keep = 0
}

# print contigs that are kept
{
    if (keep == 1)
        print $0
}
