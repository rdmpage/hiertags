# This is a perl program
# usage: 
# perl reconstr_A.pl [arguments]
# for possible arguments run the program without arguments
# more info about perl: www.perl.org


use strict;
use warnings;

my ($seed, $object, $label, $label2, $i, $j, $k, $N, $n_starting_nodes, $start, $just_reached, $n_components, $node, $node2, $n, $d, $fout, $fout2, $child, $parent, $new_dagcomponent, $r, $c, $c2, $progress, $w, $common, $label3, $max, $best, $root, $top_root, $correct, $flag, $i2, $counter, $O, $anc, $level, $current, $next, $first_label, $min, $exact, $lineage, $near, $mu, $sigma, $A, $all, $parent2, $inverted, $differentbranches, $max_h, $global_root, $loop_found, $forbidden, $use_origDAG, $use_annots, $annots_file, $use_coappearances, $coapps_file, $objn_file, $dag_file, $comment, $threshold, $help_msg);
my (@temp, @nodes, @news, @currents, @fixed_starting_nodes, @starting_nodes, @optional_starting_nodes, @sorted_neighbours, @print_order, @roots, @sorted);
my (%labels, %objects, %ll_nb, %nb, %child, %parent, %disconnected_labels, %centrality, %component, %component_sizes, %already_reached, %distance_sum, %earlier_starting_nodes, %starting_nodes_in_components, %roots, %dagcomponent, %labels_in_dagcomponent, %strength, %rank, %similarity, %objectnumber, %root, %temp, %correct, %incorrect, %ss, %neighbours, %reduced_directed_ll_nb, %suspected_ancestor, %suspected_ancestor_temp, %dontstart, %starts, %set, %nexts, %visited, %edgevalues, %maxscore_parent, %higher_level_edges, %temp_descendant, %temp_parent, %popularity, %score, %ancestors, %siblings, %zscore, %appearance, %coappeared, %zscore2, %reconstr_parent, %reconstr_child, %h, %nonglobal_roots, %suggested_parent, %suggested_component, %looped, %touched);


$threshold = 0.4;       # default value
$comment = '#';         # comment line indicator
$help_msg = "\nusage:\n\nperl reconstr_A.pl [arguments]\n\nmandatory arguments:\n\t-tags OBJECTS_WITH_TAGS_FILE (mutually exclusive with -net)\n\t-net  WEIGHTED_TAG_COAPPEARANCES_FILE (mutually exclusive with -tags)\n\t-objn OBJECTNUMBERS_FILE (only when using -net)\n\noptional arguments:\n\t-dag ORIGINAL_DAG_FILE\n\t-t LOCAL LINK WEIGHT THRESHOLD (max link weight taken into account, compared to the highest local weight, should be between 0-1. Default is 0.4)\n\t-o   OUTPUT_FILE\n\ncomment lines in the input files should begin with # (whitespaces before the # are allowed)\n\nexamples:\nperl reconstr_A.pl -tags tagged_objects_file.txt\nperl reconstr_A.pl -net tag_coappearances.txt -objn objectnumbers.txt -o output.txt -t 0.55 -dag dag_file.txt\n\n";

if (not @ARGV or @ARGV == 1 and ($ARGV[0] eq '-h' or $ARGV[0] eq '-help' or $ARGV[0] eq '--help')){
	die $help_msg;
}
for ($i = 0; $i < @ARGV; $i += 2){
	die $help_msg unless $ARGV[$i] =~ /^-/;
#	die "$ARGV[$i] requires a switch (e.g. -o)!\n" unless $ARGV[$i] =~ /^-/;
	die $help_msg if $i + 1 >= @ARGV or $ARGV[$i+1] =~ /^-/;
#	die "missing input after switch $ARGV[$i]!\n" if $i + 1 >= @ARGV or $ARGV[$i+1] =~ /^-/;
	$_ = $ARGV[$i + 1];
	die $help_msg unless $ARGV[$i] =~ /^-/;
#	die "missing switch before $_!\n" unless $ARGV[$i] =~ /^-/;
	if ($ARGV[$i] eq '-tags'){
		$use_annots = 1;
		$annots_file = $ARGV[$i + 1];
		die $help_msg if $use_coappearances;
#		die "use either -tags or -net!\n" if $use_coappearances;
	}
	elsif ($ARGV[$i] eq '-net'){
		$use_coappearances = 1;
		$coapps_file = $ARGV[$i + 1];
		die $help_msg if $use_annots;
#		die "use either -tags or -net!\n" if $use_annots;
	}
	elsif ($ARGV[$i] eq '-objn'){
		$objn_file = $ARGV[$i + 1];
		die $help_msg unless $use_coappearances;
#		die "objectnumbers are unnecessary when using -tags. Please use the -net switch with a coappearance graph input or exclude $objn_file!\n" unless $use_coappearances;
                open IN, $objn_file or die "unable to open $objn_file: $!\n";
                ($j, $k, $flag) = (0, 0, 0);
                while (<IN>){
         		if (/^\s*$comment/ or /^\s*\n/){
          		}
          		elsif (/^(\d+)\n/){
        			if ( $j != 0 and $k == 0){
                                        print STDERR "in $objn_file, first line (excluding comments and empty lines) should be a single positive integer number, the total number of objects!\n";
                                        $flag = 1;
                                }
                                unless ($k == 0){
                                        print STDERR "in $objn_file, a single number is expected only in the first (nonempty, noncomment) line!\n";
                                        $flag = 1;
                                }
                                $k++;
          		}
        		elsif (/^(\S+)\s(\d+)\n/){
         			$j++;
          		}
          		else {
         			print STDERR "invalid line in $objn_file:\n$_";
          			$flag = 1;
          		}
                }
                close IN;
                if ($j == 0){
                        print STDERR "in $objn_file, the number of objects of each tag is expected!\nformat of a line: 'tag number_of_objects'\n";
                        $flag = 1;
                }
                die if $flag;
	}
        elsif ($ARGV[$i] eq '-t'){
                $threshold = $ARGV[$i + 1];
		die $help_msg unless $threshold =~ /^\d+(\.\d+)?$/ and $threshold > -0.0000000000000005 and $threshold < 1.0000000000000005;
# 		die "threshold should be a real number between 0 and 1!\n" unless $threshold =~ /^\d+(\.\d+)?$/ and $threshold > -0.0000000000000005 and $threshold < 1.0000000000000005;
        }
	elsif ($ARGV[$i] eq '-dag'){
		$use_origDAG = 1;
		$dag_file = $ARGV[$i + 1];
	}
	elsif ($ARGV[$i] eq '-o'){
		$fout = $ARGV[$i + 1];
	}
	else {
		die $help_msg;
#		die "unknown input parameter: $ARGV[$i]!\ntype perl reconstr_zscore_A.pl for help!\n";
	}
}
die $help_msg unless $use_annots or $use_coappearances;
#die "please use either -tags or -net!\n" unless $use_annots or $use_coappearances;
die $help_msg if $use_coappearances and not $objn_file;
#die "please give an objectnumbers file when using -net!\n" if $use_coappearances and not $objn_file;


unless ($fout){
	if ($use_annots){
		$fout = $annots_file;
	}
	else {
		$fout = $coapps_file;
	}
        if ($fout =~ /\./){
                $fout =~ s/\.[^\.]+$/_tag_hierarchy_A.txt/;
        }
        else {
                $fout .= "_tag_hierarchy_A.txt";
        }
}
$fout2 = $fout;
$fout2 =~ s/tag_hierarchy_A\.txt/statistics_A.txt/;
if (-e $fout){
	$i = 2;
	do {
		if ($fout =~ /A.txt$/){
			$fout  =~ s/A\.txt$/A_$i.txt/;
			$fout2 =~ s/A\.txt$/A_$i.txt/;
		}
		else {
			$fout  =~ s/A_\d+\.txt$/A_$i.txt/;
			$fout2 =~ s/A_\d+\.txt$/A_$i.txt/;
		}
		$i++;
	} while (-e $fout);
}
open OUT2, ">".$fout2 or die "unable to open $fout2: $!\n";


print STDERR "\ninput:\n";
print OUT2   "input:\n";
print STDERR "tagged objects: $annots_file\n" if $use_annots;
print OUT2   "tagged objects: $annots_file\n" if $use_annots;
print STDERR "weighted coappearance graph: $coapps_file\nobjectnumbers: $objn_file\n" if $use_coappearances;
print OUT2   "weighted coappearance graph: $coapps_file\nobjectnumbers: $objn_file\n" if $use_coappearances;
print STDERR "link weight threshold: $threshold\n";
print OUT2   "link weight threshold: $threshold\n";
print STDERR "original DAG: $dag_file\n" if $use_origDAG;
print OUT2   "original DAG: $dag_file\n" if $use_origDAG;
print STDERR "\n";
print OUT2   "\n";



	# input

if ($use_origDAG){
	open IN, $dag_file or die "unable to open $dag_file: $!\n";
	$flag = 0;
	while(<IN>){
		if (/^\s*$comment/ or /^\s*\n/){
		}
		elsif (/^(\S+)\s(\S+)\n/){	# DAG in "child parent" edgelist format
			$parent{$1}{$2} = undef;
			$child{$2}{$1} = undef;
		
			$nb{$1}{$2} = undef;
			$nb{$2}{$1} = undef;
		}
		else {
			print STDERR "invalid line in $dag_file:\n$_";
			$flag = 1;
		}
	}
	close IN;
	die if $flag;
}


if ($use_coappearances){
	open IN, "$coapps_file" or die "unable to open $coapps_file: $!\n";
	$flag = 0;
	while (<IN>){
		if (/^\s*$comment/ or /^\s*\n/){
		}
		elsif (/^(\S+)\s(\S+)\s(\d+)\n/){
			$ll_nb{$1}{$2} = $3;
			$ll_nb{$2}{$1} = $3;
		}
		else {
			print STDERR "invalid line in $coapps_file:\n$_";
			$flag = 1;
		}
	}
	close IN;
	die if $flag;


	open IN, "$objn_file" or die "unable to open $objn_file: $!\n";
	$flag = 0;
	while (<IN>){
		if (/^\s*$comment/ or /^\s*\n/){
		}
		elsif (/^(\S+)\s(\d+)\n/){
			$objectnumber{$1} = $2;
		}
		elsif (/^(\d+)\n/){
			$O = $1;
		}
		else {
			print STDERR "invalid line in $objn_file:\n$_";
			$flag = 1;
		}
	}
	close IN;
	die if $flag;
        
        $flag = 0;
        for $label (keys %ll_nb){
                  $flag = 1 unless exists $objectnumber{$label};
        }
        die "some tags in $coapps_file do not appear in $objn_file!\n" if $flag;

        $flag = 0;
        for $label (keys %objectnumber){
                  $flag = 1 unless exists $ll_nb{$label} or not $objectnumber{$label};
        }
        die "some tags in $objn_file do not appear in $coapps_file!\n" if $flag;
}
else {
	open IN, $annots_file or die "unable to open $annots_file: $!\n";
	$flag = 0;
	while (<IN>){
		if (/^\s*$comment/ or /^\s*\n/){
		}
		elsif (/^\S+\s+(\S.*)\n/){
			@temp = split /\s/, $1;
			chomp @temp;
			undef %temp;
			for $label (@temp){
				$temp{$label} = undef;
			}
			@temp = keys %temp;
			if (@temp >= 2){
				$O += 1;
				for $label (@temp){
					$objectnumber{$label} += 1;
					for ($i = 0; $i < @temp; $i++){
						$ll_nb{$label}{$temp[$i]} += 1 unless $temp[$i] eq $label;
					}
				}
			}
		}
		else {
			print STDERR "invalid line in $annots_file:\n$_";
			$flag = 1;
		}
	}
	close IN;
	die if $flag;
}

$_ = keys %ll_nb;
print STDERR "$_ tags found\n";
print OUT2   "$_ tags found\n";
print STDERR "$O object found (having min. 2 different tags)\n";
print OUT2   "$O object found (having min. 2 different tags)\n";


for $label (keys %ll_nb){
	for $label2 (keys %{$ll_nb{$label}}){
		$strength{$label} += $ll_nb{$label}{$label2};
	}
}	




	# getting the real ancestors & siblings

if ($use_origDAG){
	for $label (keys %parent){
		undef %nexts;
		%nexts = %{$parent{$label}};
		$i = 0;
		do {
			@currents = keys %nexts;
			undef %nexts;
			$i++;
			for $current (@currents){
				$ancestors{$label}{$current} = $i unless exists $ancestors{$label}{$current};
				if (exists $parent{$current}){
					for $next (keys %{$parent{$current}}){
						$nexts{$next} = undef;
					}
				}
			}
		} while (%nexts);
	}
	for $label (keys %parent){
		for $label2 (keys %{$parent{$label}}){
			for $label3 (keys %{$child{$label2}}){
				$siblings{$label}{$label3} = undef unless $label3 eq $label;
			}
		}
	}
	for $root (keys %child){
		unless (exists $parent{$root} and $parent{$root} and keys %{$parent{$root}}){
			undef %nexts;
			$nexts{$root} = undef;
			do {
				@currents = keys %nexts;
				undef %nexts;
				for $current (@currents){
					$root{$current}{$root} = undef;;
					if (exists $child{$current}){
						for $next (keys %{$child{$current}}){
							$nexts{$next} = undef;
						}
					}
				}
			} while (%nexts);
		}
	}
}




	# finding frequent coappearing tags

for $label (keys %ll_nb){
	%neighbours = %{$ll_nb{$label}};
	@sorted_neighbours = sort weightsort keys %neighbours;
	$max = $ll_nb{$label}{$sorted_neighbours[0]};
	for $label2 (keys %{$ll_nb{$label}}){
		if ($ll_nb{$label}{$label2} >= $threshold * $max){
			$reduced_directed_ll_nb{$label}{$label2} = $ll_nb{$label}{$label2};
		}
	}
}



for $label (keys %ll_nb){
	for $label2 (keys %{$ll_nb{$label}}){
		$mu = $objectnumber{$label} * $objectnumber{$label2} / $O;
		$sigma = sqrt ($mu * ($O - $objectnumber{$label2}) / $O * ($O - $objectnumber{$label}) / ($O - 1));
		if ($sigma){
			$zscore{$label}{$label2} = ($ll_nb{$label}{$label2} - $mu) / $sigma if exists $ll_nb{$label}{$label2};
		}
		else {	# so $label or $label2 appears ALL objects
			$zscore{$label}{$label2} = 0;
		}
	}
}




	#############################
	# reconstruction - 1st part #
	#############################

open OUT, ">".$fout or die "unable to open $fout: $!\n";
print OUT "# result of the following command line:\n# reconst_A.pl @ARGV\n#\n# directed graf of tags in edge list format\n# two columns separated by space:\n# target_node source_node\n\n";
$i = 0;
($all, $exact, $lineage, $near, $inverted, $differentbranches) = (0, 0, 0, 0, 0, 0);
for $label (sort {$b cmp $a} keys %reduced_directed_ll_nb){
	$i++;
	%_ = %{$zscore{$label}};
	@_ = sort anysort keys %_;
	$parent = undef;
	for $label2 (@_){
		if (not $parent and exists $reduced_directed_ll_nb{$label}{$label2}){
			unless (exists $reduced_directed_ll_nb{$label2} and exists $reduced_directed_ll_nb{$label2}{$label}){
				$parent = $label2;
			}
		}
	}
	if ($parent){
		$reconstr_parent{$label}{$parent} = undef;
		$reconstr_child{$parent}{$label} = undef;
		$all += 1;
		if ($use_origDAG){
			if (exists $parent{$label} and exists $parent{$label}{$parent}){
				$exact += 1;
			}
			elsif (exists $ancestors{$label} and exists $ancestors{$label}{$parent}){
				$lineage += 1;
			}
			elsif (exists $siblings{$label} and exists $siblings{$label}{$parent}){
				$near += 1;
			}
			elsif (exists $ancestors{$parent} and exists $ancestors{$parent}{$label}){
				$inverted += 1;
			}
			elsif ( (not exists $ancestors{$label} or not exists $ancestors{$label}{$parent}) and (not exists $ancestors{$parent} or not exists $ancestors{$parent}{$label}) ){
				$differentbranches += 1;
			}
		}
		print OUT "$label $parent\n";
	}
	elsif ($use_origDAG and (not exists $parent{$label} or not (keys %{$parent{$label}}))){
		$exact += 1;
	}
}




	##########################
	# merging the components #
	##########################


	# finding the (reconstructed) root of each node

$max_h = 0;
for $label (keys %ll_nb){
	for $label2 (keys %{$ll_nb{$label}}){
		$h{$label} -= $ll_nb{$label}{$label2} / $strength{$label} * log ($ll_nb{$label}{$label2} / $strength{$label});
	}
	unless (exists $reconstr_parent{$label} and $reconstr_parent{$label} and keys %{$reconstr_parent{$label}}){
		$roots{$label} = undef;
		if ($max_h < $h{$label} - 0.0000000000000005){
			$max_h = $h{$label};
			$global_root = $label;
		}
	}

}
for $root (keys %roots){
	undef %nexts;
	$nexts{$root} = undef;
	do {
		@currents = keys %nexts;
		undef %nexts;
		for $current (@currents){
			$root{$current} = $root;
			if (exists $reconstr_child{$current}){
				for $next (keys %{$reconstr_child{$current}}){
					$nexts{$next} = undef;
				}
			}
		}
	} while (%nexts);
}
$_ = keys %roots;
print STDERR "number of local hierarchies is $_\n";
print OUT2   "number of local hierarchies is $_\n";


	# find a parent for the nonglobal roots

%nonglobal_roots = %roots;
delete $nonglobal_roots{$global_root};
for $root (keys %nonglobal_roots){
	%_ = %{$ll_nb{$root}};
	@_ = sort anysort keys %_;
	undef $parent;
	for $label (@_){
		unless ($parent or $root{$label} eq $root){
			$parent = $label;
			$suggested_parent{$root} = $parent;
			$suggested_component{$root} = $root{$parent};
		}
	}
}


	# detect future loops - loop handling: currently the root having the highest entropy has to find a new parent in another component

undef %looped;
for $label (keys %suggested_component){
	unless (exists $looped{$label}){
		$loop_found = 0;
		undef %touched;
		undef %nexts;
		$nexts{$label} = undef;
		do {
			@currents = keys %nexts;
			undef %nexts;
			for $current (@currents){
				if (exists $touched{$current}){
					$loop_found = 1;
				}
				else {
					$touched{$current} = undef;
					if (exists $suggested_component{$current}){
						$nexts{$suggested_component{$current}} = undef;
					}
				}
			}
		} while (%nexts);
		if ($loop_found){
			for $label2 (keys %touched){
				$looped{$label2} = $label;
				delete $suggested_parent{$label2};
				delete $suggested_component{$label2};
			}
		}
	}
}

@sorted = sort hsort keys %looped;
print STDERR "joining local hierarchies:\n";
print OUT2   "joining local hierarchies:\n";
for $root (@sorted){
	%_ = %{$ll_nb{$root}};
	@_ = sort anysort keys %_;
	undef $parent;
	for $label2 (@_){
		$forbidden = 0;
		$next = $root{$label2};
		do {
			$forbidden = 1 if $next eq $root;
			if (exists $suggested_component{$next}){
				$next = $suggested_component{$next};
			}
			else {
				$next = undef;
			}
		} while ($next and not $forbidden);
		unless ($parent or $forbidden){
			$parent = $label2;
			$suggested_parent{$root} = $parent;
			$suggested_component{$root} = $root{$parent};
		}
	}
	unless ($parent){		# if nothing else works
		$parent = $global_root;
		$suggested_parent{$root} = $parent;
		$suggested_component{$root} = $root{$parent};
	}
	print STDERR "$root $parent\n";
	print OUT2   "$root $parent\n";
}
$_ = keys %suggested_parent;
print STDERR "\nnumber of new links is $_\n";
print OUT2   "\nnumber of new links is $_\n";




for $root (sort {$b cmp $a} keys %suggested_parent){
	$all += 1;
	if ($use_origDAG){
		if (exists $parent{$root} and exists $parent{$root}{$suggested_parent{$root}}){
			$exact += 1;
		}
		elsif (exists $ancestors{$root} and exists $ancestors{$root}{$suggested_parent{$root}}){
			$lineage += 1;
		}
		elsif (exists $siblings{$root} and exists $siblings{$root}{$suggested_parent{$root}}){
			$near += 1;
		}
		elsif (exists $ancestors{$suggested_parent{$root}} and exists $ancestors{$suggested_parent{$root}}{$root}){
			$inverted += 1;
		}
		elsif ( (not exists $ancestors{$root} or not exists $ancestors{$root}{$suggested_parent{$root}}) and (not exists $ancestors{$suggested_parent{$root}} or not exists $ancestors{$suggested_parent{$root}}{$root}) ){
			$differentbranches += 1;
		}
	}
	print OUT "$root $suggested_parent{$root}\n";
}
close OUT;
if ($use_origDAG){
	$_ = $exact + $lineage + $near + $inverted;
	print STDERR "exact = $exact\nancestor = $lineage\ninverted = $inverted\nsibling = $near\ntotal related = $_\ndiffbranches (incl. siblings) = $differentbranches\n";
	print OUT2   "exact = $exact\nancestor = $lineage\ninverted = $inverted\nsibling = $near\ntotal related = $_\ndiffbranches (incl. siblings) = $differentbranches\n";
}
$_ = (keys %ll_nb) - 1;
print STDERR "total number of found edges = $all (max possible $_)\n\noutput written to $fout\n\n";
print OUT2   "total number of found edges = $all (max possible $_)\n";
close OUT2;




sub weightsort {$neighbours{$b} <=> $neighbours{$a}}

sub hsort {
	if (abs ($h{$b} - $h{$a}) >= 0.0000000000000005){
		$h{$b} <=> $h{$a}
	}
	else {
		$b cmp $a;
	}
}

sub anysort {
	if (abs ($_{$b} - $_{$a}) >= 0.0000000000000005){
		$_{$b} <=> $_{$a};
	}
	else {
		$b cmp $a;
	}
}


