# This is a perl program
# usage: 
# perl reconstr_B.pl [arguments]
# for possible arguments run the program without arguments
# more info about perl: www.perl.org


use strict;
use warnings;
use POSIX;

my ($label, $label2, $label3, $i, $j, $k, $fout, $fout2, $fout3, $parent, $flag, $O, $current, $next, $exact, $lineage, $near, $mu, $sigma, $threshold, $max_w, $max_z, $max_d, $c, $all, $inverted, $differentbranches, $use_origDAG, $use_annots, $annots_file, $use_coappearances, $coapps_file, $objn_file, $dag_file, $comment, $help_msg, $i0, $component);
my (@temp, @currents, @temp2, @componentlist);
my (%ll_nb, %nb, %child, %parent, %centrality, %strength, %objectnumber, %temp, %nexts, %visited, %score, %ancestors, %siblings, %zscore, %temp2, %ll_nb_filtered, %descendantsvotes, %child_reconstr, %sum, %sum2, %component, %lambda, %lambda2);


$threshold = 10;       # default value
$comment = '#';        # comment line indicator
$help_msg = "\nusage:\n\nperl reconstr_B.pl [arguments]\n\nmandatory arguments:\n\t-tags OBJECTS_WITH_TAGS_FILE (mutually exclusive with -net)\n\t-net  WEIGHTED_TAG_COAPPEARANCES_FILE (mutually exclusive with -tags)\n\t-objn OBJECTNUMBERS_FILE (only when using -net)\n\noptional arguments:\n\t-dag ORIGINAL_DAG_FILE\n\t-z   Z-SCORE_THRESHOLD_VAUE (a nonnegative real number; default is 10)\n\t-o   OUTPUT_FILE\n\ncomment lines in the input files should begin with # (whitespaces before the # are allowed)\n\nexamples:\nperl reconstr_B.pl -tags tagged_objects_file.txt\nperl reconstr_B.pl -net tag_coappearances.txt -objn objectnumbers.txt -o output.txt -z 8.67 -dag dag_file.txt\n\n";

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
	elsif ($ARGV[$i] eq '-z'){
		$threshold = $ARGV[$i + 1];
		die $help_msg unless $threshold =~ /^\d+(\.\d+)?$/;
#		die "$threshold is not a valid z-score threshold!\n" unless $threshold =~ /^\d+(\.\d+)?$/;
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
#		die "unknown input parameter: $ARGV[$i]!\ntype perl reconstr_flickr_simple_B.pl for help!\n";
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
                $fout =~ s/\.[^\.]+$/_tag_hierarchy_B.txt/;
        }
        else {
                $fout .= "_tag_hierarchy_B.txt";
        }
}
$fout2 = $fout;
$fout2 =~ s/tag_hierarchy_B\.txt/statistics_B.txt/;
$fout3 = $fout;
$fout3 =~ s/tag_hierarchy_B\.txt/centralities_B.txt/;
if (-e $fout){
	$i = 2;
	do {
		if ($fout =~ /B.txt$/){
			$fout  =~ s/B\.txt$/B_$i.txt/;
			$fout2 =~ s/B\.txt$/B_$i.txt/;
		}
		else {
			$fout  =~ s/B_\d+\.txt$/B_$i.txt/;
			$fout2 =~ s/B_\d+\.txt$/B_$i.txt/;
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
print STDERR "z-score threshold = $threshold\n";
print OUT2   "z-score threshold = $threshold\n";
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
                  $flag = 1 unless exists $ll_nb{$label} or not $objectnumber{$label};;
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


	# centrality calculation

if (-e $fout3){
	open IN, $fout3 or die "unable to open $fout3: $!\n";
	while (<IN>){
		if (/^(\S+)\s+(\S+)\n/){
			$centrality{$1} = $2;
		}
		elsif (/^residual=(\S+)\n/){
			$c = $1;
		}
	}
	close IN;
	print STDERR "centralities read from $fout3\n";
	print OUT2   "residual = $c (centralities read from $fout3)\n";
}
else {
	for $label (keys %ll_nb){
		for $label2 (keys %{$ll_nb{$label}}){
			$ll_nb_filtered{$label}{$label2} = $ll_nb{$label}{$label2} if $zscore{$label}{$label2} >= $threshold or ($ll_nb{$label}{$label2} / $objectnumber{$label2} > 0.5 or $ll_nb{$label}{$label2} / $objectnumber{$label} > 0.5);
		}
		$strength{$label} = 0;
		for $label2 (keys %{$ll_nb{$label}}){
			$strength{$label} += $ll_nb{$label}{$label2} if $zscore{$label}{$label2} >= $threshold or ($ll_nb{$label}{$label2} / $objectnumber{$label2} > 0.5 or $ll_nb{$label}{$label2} / $objectnumber{$label} > 0.5);
		}
		$centrality{$label} = 0 unless exists $ll_nb_filtered{$label};
	}


	print STDERR "centrality calculation\n";
	@_ = keys %ll_nb_filtered;
	
		# finding connected components in %ll_nb_filtered

	$i0 = 0;
	do {
		for ($i = $i0; $i < @_ and exists $visited{$_[$i]}; $i++){
		}
		$i0 = $i;
		
		$label = $_[$i];
		undef %nexts;
		$nexts{$label} = undef;
		$visited{$label} = undef;
		$component{$label} = $label;
		push @componentlist, $label;
		while (%nexts){
			@currents = keys %nexts;
			undef %nexts;
			for $current (@currents){
				for $next (keys %{$ll_nb_filtered{$current}}){
					unless (exists $visited{$next}){
						$nexts{$next} = undef;
						$visited{$next} = undef;
						$component{$next} = $label;
					}
				}
			}
		}
	} while ((keys %visited) < (keys %ll_nb_filtered));

		# calculating the first eigenvector by power iteration. Normalisation is done separately for each connected component.

	for $component (@componentlist){
		$sum{$component} = 0;
	}
	undef %temp;
	for $label (@_){
		$temp{$label} = $strength{$label};
		$sum{$component{$label}} += abs $strength{$label};
	}
	for ($i = 0; $i < 100; $i++){	# IMPORTANT:	we will assume that there is an EVEN number of iterations (because we potentially need to deal with an oscillating case, for which it is more convenient and the code assumes an even number)
		unless ( (POSIX::floor ($i / 2)) * 2 < $i - 0.0000000000000005 ){	# $i is even
			for $component (@componentlist){
				$sum2{$component} = 0;
			}
			for $label (@_){
				$temp2{$label} = 0;
				for $label2 (keys %{$ll_nb_filtered{$label}}){
					$temp2{$label} += $temp{$label2} * $ll_nb_filtered{$label}{$label2} / $sum{$component{$label}};
				}
				$sum2{$component{$label}} += abs $temp2{$label};
			}
		}
		else {									# $i is odd
			for $component (@componentlist){
				$sum{$component} = 0;
			}
			for $label (@_){
				$temp{$label} = 0;
				for $label2 (keys %{$ll_nb_filtered{$label}}){
					$temp{$label} += $temp2{$label2} * $ll_nb_filtered{$label}{$label2};	# because we norm the iterated vector only after every second iteration, the extra length gained in two iterations is always lambda_max^2 (assuming smaller eigenvalues have already died out), while in the oscillating case one iteration increases the length by a factor of a mixture containing the components of the two affected eigenvectors etc.
				}
				$sum{$component{$label}} += abs $temp{$label};
			}
		}
		print STDERR "$i/100 ";
	}
	$c = 0;
	for $label (@_){
		$centrality{$label} = $temp{$label} / $sum{$component{$label}};
		$c += abs ($centrality{$label} - $temp2{$label} / $sum2{$component{$label}});
	}
	print STDERR "\nresidual = $c\n";


		# if the solution oscillates, i.e. there is also a -lambda_max eigenvalue: cancel out the "negative" eigenvector component
		# it can happen that (some components of) the coappearance graph is almost bipartite (e.g. a tree + 1 weak link), then there is an almost (-1)*lambda_max eigenvalue, which slows down the convergence but which is not totally cancelled out by the oscillation solving procedure. In this case, the residual (the approximated error) can be much higher than 10^-17. If you are not satisfied, apply more iterations.

	if ($c > 0.0000000001){
		for $component (@componentlist){
			$lambda{$component} = sqrt $sum{$component};		# $sum is lambda_max^2 (lambda is the largest eigenvalue), but is is more clear to use lambda-names a few lines later
			$lambda2{$component} = $sum{$component};
		}
		for $component (@componentlist){
			$sum{$component} = 0;
		}
		for $label (@_){
			$centrality{$label} = 0.5 * ($temp{$label} / $lambda2{$component{$label}} + $temp2{$label} / $lambda{$component{$label}});
			$sum{$component{$label}} += abs $centrality{$label};
		}

		for $component (@componentlist){		# check the residual
			$sum2{$component} = 0;
		}
		for $label (@_){
			$temp2{$label} = 0;
			for $label2 (keys %{$ll_nb_filtered{$label}}){
				$temp2{$label} += $centrality{$label2} * $ll_nb_filtered{$label}{$label2} / $sum{$component{$label}};
			}
			$sum2{$component{$label}} += abs $temp2{$label};
		}
		$c = 0;
		for $label (@_){
			$c += abs ($temp2{$label} / $sum2{$component{$label}} - $centrality{$label} / $sum{$component{$label}});
			$centrality{$label} /= $sum{$component{$label}};
		}
		print STDERR "\npotentially oscillating solution fixed\n";
		print STDERR "residual = $c\n";
	}
	open OUT3, ">".$fout3 or die "unable to open $fout3: $!\n";
	for $label (reverse sort centralitysort @_){
		print OUT3 "$label\t$centrality{$label}\n";
	}
	print OUT2   "residual = $c\n";
	print OUT3   "residual=$c\n";
	close OUT3;

	print STDERR "done\n";
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
}


	# reconstruction

($all, $exact, $lineage, $near, $inverted, $differentbranches) = (0, 0, 0, 0, 0, 0);
open OUT, ">".$fout or die "unable to open $fout: $!\n";
print OUT "# result of the following command line:\n# reconst_B.pl @ARGV\n#\n# directed graf of tags in edge list format\n# two columns separated by space:\n# target_node source_node\n\n";
for $label (sort centralitysort keys %ll_nb){
	undef %descendantsvotes;
	undef %nexts;
	$nexts{$label} = undef;
	do {
		@currents = keys %nexts;
		undef %nexts;
		for $current (@currents){
			for $next (keys %{$child_reconstr{$current}}){
				for $label2 (keys %{$ll_nb{$next}}){
					$descendantsvotes{$label2} += $zscore{$label2}{$next} if exists $ll_nb{$label}{$label2} and $centrality{$label2} > $centrality{$label} and ($zscore{$label}{$label2} >= $threshold or $ll_nb{$label}{$label2} / $objectnumber{$label} > 0.5) and ($zscore{$next}{$label2} >= $threshold or $ll_nb{$next}{$label2} / $objectnumber{$next} > 0.5);
				}
				$nexts{$next} = undef; 
			}
		}
	} while (%nexts);
	
	undef %score;
	for $label2 (keys %{$ll_nb{$label}}){
		if ($centrality{$label2} > $centrality{$label} and ($zscore{$label}{$label2} >= $threshold or $ll_nb{$label}{$label2} / $objectnumber{$label} > 0.5)){
			$descendantsvotes{$label2} += $zscore{$label}{$label2};
			$score{$label2} = $descendantsvotes{$label2};
		}
	}

	$parent = undef;
	if (%score){
		$parent = (sort scoresort keys %score)[0];
	}
	if ($parent){
		$child_reconstr{$parent}{$label} = undef;
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
close OUT;
if ($use_origDAG){
	$_ = $exact + $lineage + $near + $inverted;
	print STDERR "exact = $exact\nancestor = $lineage\ninverted = $inverted\nsibling = $near\ntotal related = $_\ndiffbranches (incl. siblings) = $differentbranches\n";
	print OUT2   "exact = $exact\nancestor = $lineage\ninverted = $inverted\nsibling = $near\ntotal related = $_\ndiffbranches (incl. siblings) = $differentbranches\n";
}
$_ = (keys %ll_nb) - 1;
print STDERR "\ntotal number of found edges = $all (max possible $_)\n\noutput written to $fout\n\n";
print OUT2   "\ntotal number of found edges = $all (max possible $_)\n";
close OUT2;



sub centralitysort {
	if (abs ($centrality{$a} - $centrality{$b}) >= 0.0000000000000005){
		$centrality{$a} <=> $centrality{$b};
	}
	else {
		$b cmp $a;
	}
}


sub scoresort {
	if (abs ($score{$b} - $score{$a}) >= 0.0000000000000005){
		$score{$b} <=> $score{$a};
	}
	else {
		$b cmp $a;
	}
}
