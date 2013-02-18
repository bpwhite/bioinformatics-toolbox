# A simulation of pseudo-biological entities that can be selected
# for certain parameters, including but no limited to: information,
# structure, growth, and replication.

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

use SimulationObject;
use Point;
use Force;
use Grid;
use Simulation;

use strict;
use warnings;


SimulationObject->logger_prefix('simulation_log');
SimulationObject->clear_log_file;

# Alphabet parameters
Alphabet->mean_alphabet_size(50);

# Letter parameters
Letter->letter_distribution('gaussian');
Letter->mean_letter_radius(25);
Letter->mean_letter_num_centroids(5);
Letter->letter_std_dev_rate(0.50);

# Centroid parameters
Centroid->centroid_distribution('gaussian');
Centroid->centroid_std_dev_rate(0.05);
Centroid->mean_steric_radius(0.25);
Centroid->mean_electrostatic_radius(0.25);
Centroid->mean_hydrogen_radius(0.25);
Centroid->mean_vanderwaals_radius(0.25);
Centroid->mean_covalent_radius(0.25);


my $simulation1 = Simulation->new( grid_xmax => 100, grid_ymax => 100);
