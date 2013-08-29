#!/usr/bin/perl
# Insert a large CSV file into a remote sql database.
#
# Copyright (c) 2013, Bryan White, bpcwhite@gmail.com

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

#!/usr/bin/perl

# DBI is the standard database interface for Perl
# DBD is the Perl module that we use to connect to the MySQL database
use FindBin;
use lib "$FindBin::Bin/libs/";
use General::Arguments;


use DBI;

use strict;
use warnings;

my $params = General::Arguments->new(	arguments_v => \@ARGV,
										option_defs => {'-database' 	=> '', 				# Database name
														'-location'		=> 'localhost',		# IP address/hostname of database
														'-port'			=> '3306',			# DB port
														'-user'			=> '',				# Username
														'-password'		=> '',				# User password	
														'-table'		=> '',				# Table to input data to
														'-csv-file'		=> '',				# CSV file to upload
													}
													);

my $database 	=  $params->options->{'-database'};
my $location 	=  $params->options->{'-location'};
my $port		=  $params->options->{'-port'};
my $user 		=  $params->options->{'-user'};
my $password 	=  $params->options->{'-password'};
my $table 		=  $params->options->{'-table'};
my $csv_file	=  $params->options->{'-csv-file'};

print "Opening $csv_file \n";
open CSV, "< $csv_file";
my @csv_file = <CSV>;
close CSV or die "Couldn't close $csv_file!\n";

my $dsn = "DBI:mysql:database=$database;host=$location;port=$port";
my $dbh = DBI->connect($dsn, $user, $password);

my $line_counter = 1;
foreach my $line (@csv_file) {
	if ($line_counter == 1) {
		$line_counter++;
		next;
	}
	# print $line."\n";
	# $line =~ s/"//g; # Strip quotes
	# print $line."\n";
	my @records = split(/",/,$line);
	my $num_records = scalar @records;
	my $insert_string = '';
	my $record_counter = 0;
	foreach my $split (@records) {
		$record_counter++;
		$split =~ s/"//g; # Strip double quotes
		$split =~ s/'//g; # Strip single quotes
		
		$split = 'NA' if $split eq ''; # Fill blanks with NA
		# print $split."\n";
		if ($record_counter != $num_records) {
			$insert_string .= "'".$split."',";
		} else {
			$insert_string .= $split;
		}
		
	}
	$insert_string =~ s/,+$//;
	
	# my $sth = $dbh->prepare("SELECT * FROM foo WHERE bla");
	# foreach my 
	my $query = "INSERT into $table values  ( NULL, ".$insert_string." )";

	my $statement = $dbh->prepare($query);
	# eval {
	$statement->execute() or die "$@\n";
	# }
	# if($@) {
		# die "Could not insert record $@\n";
	# }
	$line_counter++;
}

