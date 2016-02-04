# Merge genbank formatted files downloaded using the genbank converter
# Created by Bryan White, 2016
#
#



open (TAXALIST, '<'.$taxa_file);
@taxa_list = <TAXALIST>;
close(TAXALIST);
