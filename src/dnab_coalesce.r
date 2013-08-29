# Delimits species using the SPLITS package in R.

# Copyright (c) 2013, Bryan White

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

require(methods,quietly=TRUE)
require(subplex,quietly=TRUE)
require(ape,quietly=TRUE)
require(splits,quietly=TRUE)
require(geiger,quietly=TRUE)
source("spec.list.1.1.r")

current_tree <- read.nexus("your_tree.tre", tree.names = NULL)
current_tree <- multi2di(current_tree, random = FALSE)
gmyc_results <- gmyc(current_tree, method="single", interval=c(0, 10))
print(summary(gmyc_results))

species_list = spec.list(gmyc_results)
write.table(species_list,file="output.txt",row.names = FALSE, sep=",")
tip_colors = c()
for(i in 1:length(current_tree$tip.label)) {
	for(j in 1:length(species_list$GMYC_spec)) {
		if(current_tree$tip.label[i] == toString(species_list$sample_name[j])) {
			tip_colors = c(tip_colors,species_list$GMYC_spec[j])
			current_tree$tip.label[i] = paste("[",species_list$GMYC_spec[j],"]_",current_tree$tip.label[i],sep="")
			break
		}
	}
}
print_to_pdf = TRUE
if(print_to_pdf) {
	pdffilename = "output_tree.pdf"
	pdf(pdffilename,
		width = 10, height = 7.5, 
		pointsize = 12, bg = "white")
}

	plot.phylo(	current_tree, "p", font=1, cex=0.05,
				show.tip.label=TRUE,
				tip.color = tip_colors)

axisPhylo()
dev.off()
