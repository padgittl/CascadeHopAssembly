#!/usr/bin/perl -w
use Graph;
use Graph::Directed;
use GraphViz;
use strict;
$|=1;

my $usage = "Usage:\n$0 <lastz results> <contig length file> [parsed purge haplotigs file]\n";

my $minCoverage = 20;
my $minScoreDensity = 20;
my $minContinuity = 20;
my $minIdentity = 20;
my $maxBigPerc = 40;
my $contigLastzResultsFile = $ARGV[0] or die $usage;
my $contigLengthFile = $ARGV[1] or die $usage;
my $purgeHaplotigsFile = $ARGV[2] ? $ARGV[2] : "";
my $contigLength = readContigLengthFile($contigLengthFile);
processLASTZContigGraph($contigLastzResultsFile,$contigLength,$purgeHaplotigsFile);

sub readContigLengthFile {
    my($contigLengthFile) = @_;
    my %contigLength;
    open(FILE,$contigLengthFile) or die "Could not open $contigLengthFile\n";
    while(<FILE>) {
	chomp;
	my($contig,$length) = split();
	$contigLength{$contig} = $length;
    }
    return \%contigLength;
}

sub processLASTZContigGraph {
    my($contigLastzResultsFile,$contigLength,$purgeHaplotigsFile) = @_;
    my $G = Graph::Undirected->new;
    my $DG = Graph::Directed->new;
    open(FILE,$contigLastzResultsFile) or die "could not open $contigLastzResultsFile\n";
    
    while(<FILE>) {
	chomp;
	# assumed that contig2 is larger
	# start,end are positions in contig2
	my($contig1,$contig2,$start,$end,$strand,$score,$scoreDensity,$coverage,$percentIdentity,$percentContinuity,$bigPerc) = split();
	if(($coverage >= $minCoverage)&&($scoreDensity >= $minScoreDensity)&&($percentContinuity >= $minContinuity)&&($percentIdentity >= $minIdentity)) {
	    if($bigPerc <= $maxBigPerc) {
		$G->add_edge($contig1,$contig2);
		$DG->add_edge($contig1,$contig2);
	    } else {
		#print "lost $contig1 $contig2 $scoreDensity vs $bigPerc\n";
	    }
	}
    }
    close(FILE);
    if($purgeHaplotigsFile) {
	open(PHF,$purgeHaplotigsFile) or die "Could not open $purgeHaplotigsFile\n";
	while(<PHF>) {
	    chomp;
	    my($contig1,$contig2,$start,$end,$strand,$score,$scoreDensity,$coverage,$percentIdentity,$percentContinuity,$bigPerc) = split();
	    if(($coverage >= $minCoverage)&&($scoreDensity >= $minScoreDensity)&&($percentContinuity >= $minContinuity)&&($percentIdentity >= $minIdentity)) {
		if($bigPerc <= $maxBigPerc) {
		    $G->add_edge($contig1,$contig2);
		    $DG->add_edge($contig1,$contig2);
		} else {
		    #print "lost $contig1 $contig2 $scoreDensity vs $bigPerc\n";
		}
	    }
	}
	close(PHF)
    }
    my @CC = sort{scalar(@{$b}) <=> scalar(@{$a})} $G->connected_components;
    my %map;
    my %isAssociate;
    # for each connected component
    for my $cc (@CC) {
	# get vertex/contig with largest IN DEGREE
	my($vSet,$d) = largestInDegree($cc,$DG,$contigLength);
	for my $v (@{$vSet}) {
	    for my $w ($DG->neighbors($v)) {
		# for all w->v
		if($DG->has_edge($w,$v)) {		    
		    # directly connected to a large primary contig.
		    #print "$v<-$w\tAssociate\n";		    
		    unless($isAssociate{$v}) {
			push(@{$map{$v}},$w);
			$isAssociate{$w} = $v;
		    }
		} elsif($DG->has_edge($v,$w)) {
		    die "EXCEPTION FOUND: Found edge going out of a max in-degree node. $v -> $w\n";
		} else {
		    die "EXEPTION: UNKNOWN.\n";
		}
		for my $x ($DG->neighbors($w)) {
		    # for all x->w->v
		    if($x ne $v) {
			if($DG->has_edge($x,$w)) {
			    # SECONDARY CONTIG
			    #print "$v<-$w<-$x\tSecondary\n";
			} elsif($DG->has_edge($w,$x)) {
			    # SECONDARY CONTIG
			    #print "$v<-$w->$x\tSecondary\n";			    
			} else {
			    die "EXCEPTON(3) $v $w $x\n";
			}				
			for my $y ($DG->neighbors($x)) {
			    # for all y->x->w->v
			    if($y ne $w) {
				if($DG->has_edge($y,$x)) {
				    #print "$v<-$w--$x<-$y\tSecondary-associate\n";
				    # it should be that X is the primary, and Y is the associate
				    unless($isAssociate{$x}) {
					push(@{$map{$x}},$y);
					$isAssociate{$y} = $x;
				    }
				} else {
				    # it should be that Y is the primary and X is the associate
				    #die "EXCEPTION(4): $v<-$w--$x->$y\n";
				    unless($isAssociate{$y}) {
					push(@{$map{$y}},$x);
					$isAssociate{$x} = $y;
				    }
				}		    
				for my $z ($DG->neighbors($y)) {
				    # for all z->y->x->w->v
				    if($z ne $x) {
					#die "EXCEPTION(5): $v $w $x $y $z\n";
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    my $COUNT = 0;
    my $totalAssociateLength = 0;
    my $clusterFile = "clusters_o".$maxBigPerc."_dens".$minScoreDensity.".txt";
    open(CF,">$clusterFile");
    for my $v (sort {scalar(@{$map{$b}}) <=> scalar(@{$map{$a}})} keys %map) {
	unless($isAssociate{$v}) {
	    $COUNT++;
	    # add 1 to size to include primary
	    print CF "Cluster-$COUNT\tsize=",scalar(@{$map{$v}})+1,"\n";
	    print CF "$v\tPrimary (",$contigLength->{$v},"bp)\n";
	    for my $w (@{$map{$v}}) {
		print CF "$v<-$w\tAssociate (",$contigLength->{$w},"bp)\n";
		$totalAssociateLength += $contigLength->{$w};
	    }
	    print CF "\n";
	}
    }
    print CF "##Total Primary: $COUNT\n";
    print CF "##Associate\t$totalAssociateLength bp\n";
}

sub largestInDegree {
    # takes a connected component, and a graph object
   my($cc,$DG,$contigLengths) = @_;
   my $maxDegree = 0;
   my $maxV;
   for my $v (@{$cc}) {
       my $d = $DG->in_degree($v);
       #print ">>$v $d\n";
       if($d > $maxDegree) {
	   $maxV = $v;
	   $maxDegree = $d;
       }
   }
   my @maxDegreeList;
   for my $v (@{$cc}) {
       my $d = $DG->in_degree($v);
       if($d == $maxDegree) {
	   push(@maxDegreeList,$v);
       }
   }
   if(@maxDegreeList > 1) {
       my @candidates;
       for my $v (@maxDegreeList) {
	   if($DG->out_degree($v) == 0) {
	       push(@candidates,$v);
	   }
       }
       if(@candidates > 1) {
	   #my $maxLength = 0;
	   #my $maxW;	       
	   #for my $w (@candidates) {
	   #    if(defined($maxLength)) {
	   #	   if($contigLengths->{$w} > $maxLength) {
	   #	       $maxLength = $contigLengths->{$w};
	   #	       $maxW = $w;
	   #	   }
	   #    } else {
	   #	   die "No length found for $w\n";
	   #    }
	   #}
	   return(\@candidates,$maxDegree);
       } elsif(@candidates == 0) {
	   die "No candiates with 0 out-degree @maxDegreeList";
       }
       # must only be one candidate:
       $maxV = shift(@candidates);
       return(\@candidates,$maxDegree);
   } else {
       $maxV = shift(@maxDegreeList);
       # if max degree node has out-edges, they must be to the real primary
       if($DG->out_degree($maxV)) {
	   my @candidates;
	   # collect all nodes that are pointed to
	   for my $w ($DG->neighbors($maxV)) {
	       if($DG->has_edge($maxV,$w)) {
		   push(@candidates,$w);
	       }
	   }
	   if(@candidates == 1) {
	       my $maxW = shift(@candidates);
	       if($DG->out_degree($maxW)) {
		   die "EXCEPTION: no max found! @{$cc}\n";
	       }
	       my @vSet = ($maxW);
	       return(\@vSet,1);
	   }
       }
       my @vSet = ($maxV);
       return(\@vSet,$maxDegree);
   }
   die "Expection found in maxInDegree.\n";
}

sub goodLinkage {
    my($contig1,$contig2) = @_;
    my %filter;
    # assume $contig1 is shorter
    $filter{"006536F"}{"002388F"} = 1;
    $filter{"006617F"}{"007332F"} = 1;
    $filter{"011680F"}{"002183F"} = 1; # couldn't even run mummerplot on this one. 
    $filter{"011410F"}{"000817F"} = 1;
    $filter{"011410F"}{"010975F"} = 1;
    $filter{"011167F"}{"004666F"} = 1;
    $filter{"011167F"}{"000120F"} = 1;
    $filter{"011711F"}{"000126F"} = 1;
    $filter{"009469F"}{"001255F"} = 1;
    $filter{"005916F"}{"002747F"} = 1;
    $filter{"006512F"}{"001482F"} = 1;
    $filter{"004148F"}{"003982F"} = 1;
    $filter{"011926F"}{"002062F"} = 1;
    $filter{"010792F"}{"001739F"} = 1; # overlaps another alignment of 010792F to 006300F 
    $filter{"010702F"}{"002909F"} = 1; # no alignment to plot
    $filter{"008942F"}{"000027F"} = 1;
    $filter{"009541F"}{"004171F"} = 1; # looks repetive, overlaps alignment of 009541F to 008550F    
    $filter{"011721F"}{"001404F"} = 1;
    $filter{"002542F"}{"002540F"} = 1; # purge haplotig chooses shorter as primary contig
    $filter{"002719F"}{"004797F"} = 1; # purge haplotig chooses shorter as primary contig
    $filter{"003097F"}{"003357F"} = 1; # purge haplotig chooses shorter as primary contig
    $filter{"002632F"}{"002882F"} = 1; # purge haplotig chooses shorter as primary contig
    $filter{"003804F"}{"003808F"} = 1;  # purge haplotig chooses shorter as primary contig
    $filter{"007419F"}{"007631F"} = 1;  # purge haplotig chooses shorter as primary contig
    $filter{"002448F"}{"003353F"} = 1;  # purge haplotig chooses shorter as primary contig
    $filter{"003354F"}{"004249F"} = 1;  # purge haplotig chooses shorter as primary contig
    $filter{"004582F"}{"004620F"} = 1;  # purge haplotig chooses shorter as primary contig
    $filter{"002210F"}{"002376F"} = 1;  # purge haplotig chooses shorter as primary contig

    if($filter{$contig1}{$contig2}) {
	return 0;
    }
    return 1;
}
