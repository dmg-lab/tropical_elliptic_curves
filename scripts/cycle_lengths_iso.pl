use strict;
use warnings;
use Benchmark qw(:all);
use POSIX qw(strftime);
use FileHandle;
use Sys::Hostname;
use File::Basename qw( dirname );
use IO::Uncompress::UnXz qw(unxz $UnXzError) ;

use application 'polytope';
my $codename="collect_graphs_";

###############################################################################
###############################################################################
### cycle_lengths_iso.pl
###
### This script runs over the output of mptopcom or TOPCOM and produces the
### graphs as described in ??. Furthermore it collects representatives of these
### graphs up to graph isomorphism.
###
### Sample usage:
### polymake --script cycle_lengths_iso.pl input.dat regular_unimodular.dat.xz
### Run in parallel as follows:
### ls -1 regular_unimodular.0*.dat.xz | parallel --eta -v --progress --bar\
###    -n 1 -j5 polymake --script cycle_lengths_iso.pl input.dat {}
###
###############################################################################
###############################################################################


die "usage: $codename MPTOPCOM_INPUT MPTOPCOM_OUTPUT<n>\n" unless scalar(@ARGV)==2;

# Set some file names, read original points to later make sure that the order
# returned from polymake agrees with them.
my $input_file = $ARGV[0];
my $input = `cat $input_file`;
$input =~ s/\s//g;
my ($pts, $gp) = $input =~ m/(.*\]\])(\[\[.*)/;
$pts = new Matrix(eval($pts));
my $tdata=$ARGV[1];
my $output_dir = dirname(`realpath $tdata`);
print "Will write output to $output_dir\n";
my $root = `realpath $tdata`;
chomp $root;
$root =~ s/\.dat\.xz//;
my $output_file="$root.graphs";
my $output_classes = "$output_file.classes";
my $log_file="$output_file.log";
# Someone else is working on this file.
if(-e $output_file || -e $log_file || -e "$output_file.xz"){
   exit;
}

# mark start
my $host = hostname();
my $now = strftime "%a %b %e %H:%M:%S %Y", localtime();
my $LOG = FileHandle->new("> $log_file");
die "cannot write $log_file\n" unless defined($LOG);
$LOG->autoflush(1);
print $LOG "started $codename @ $host on $now\n";
print $LOG "reading $tdata\n";
# open my $TDATA,  "xzcat $tdata|" or die "cannot read $tdata";
my $TDATA = IO::Uncompress::UnXz->new($tdata)
    or die "IO::Uncompress::UnXz failed: $UnXzError\n";

print $LOG "writing $output_file\n";
my $OUT = FileHandle->new("> $output_file");
$OUT->autoflush(1);
die "cannot write $output_file\n" unless defined($OUT);


my $P = scale(simplex(3),2);
my $CP = cayley_polytope($P,$P);
my $cayley_points = $CP->LATTICE_POINTS;

# Check that the points are the same as given to mptopcom.
if($cayley_points->minor(All,~[5]) != $pts){
   die "Points do not agree with input points from mptopcom!";
}



###############################################################################
# Function for computing graph from a triangulation.
sub tropical_curve_from_subdivision{
   my($pts, $mc, $index_map) = @_;
   my $toblerone_cells = new Set<Set<Integer>>;

   foreach my $cell (@{$mc}){
      my $set1 = new Set<Vector<Rational>>;
      my $set2 = new Set<Vector<Rational>>;
      foreach my $vertex (@{$pts->minor($cell, All)}){
         if ($vertex->[4] == 1){$set1 += $vertex} else {$set2 += $vertex};}   
      # There are only three possible cases for the convex hulls of $set1 and
      # $set2: a point, a line, or a triangle. In any case the elements of
      # $set1 and $set2 are the vertices of a maximal cell, hence these will be
      # vertices of the convex hulls as well. We only want those cases where
      # the Minkowski sum of these convex hulls is a toblerone, i.e. a prism
      # over a triangle. This happens for the case of a line and a triangle.
      # Thus the sizes of $set1 and $set2 must be 2 and 3, or the other way
      # around. Thus we can just test that the product of the sizes is 6, this
      # already excludes all other cases. To see that this is sufficient,
      # notice that $set1 and $set2 now form a disjoint subdivision of the
      # vertices of a 4-dim simplex.
      if($set1->size()*$set2->size() == 6){
         my $sum_vert = new Set<Vector<Integer>>();
         foreach my $v1 (@$set1){
            foreach my $v2 (@$set2){
               my $s = $v1+$v2;
               $s->[0] = 1;
               $sum_vert += $s->slice(range(0,3));
            }
         }
         my $cell_indices = new Set<Integer>;
         foreach my $vertex (@{$sum_vert}){
            $cell_indices += $index_map->{$vertex};
         }
         $toblerone_cells += $cell_indices
      }
   }

   my $ordered_toblerones_cells = new Array<Set<Int>>($toblerone_cells);
   my $edges = new Set<Array<Int>>;
   for (my $i=0; $i<15; ++$i) {
      for (my $j=$i+1; $j<16; ++$j) {
         my $intersection = $ordered_toblerones_cells->[$i]*$ordered_toblerones_cells->[$j];
         if ($intersection->size() == 4) {$edges += [$i,$j];}
      }
   }
   return $edges;
}


###############################################################################
# Compute point set and compare it to mptopcom input points to make sure that
# these are in the right order.
my $minkowski = minkowski_sum($P, $P);
my $n_points = ($minkowski->LATTICE_POINTS)->rows();
my $points = $minkowski->LATTICE_POINTS;
my $ptsindices = new Map<Vector<Integer>, Int>();
for(my $i=0; $i<$points->rows(); $i++){
   $ptsindices->{$points->row($i)} = $i;
}
my @collected_graphs = ();

my %histogram = ();
my $classes = new Map<Int, GraphAdjacency>();
my $nclasses = 0;


# start the clock
my $t0=Benchmark->new();


my $id=0;
while(my $line=$TDATA->getline()) {
   chomp $line;
   $line =~ s/^.*(\{\{.*\}\}).*/$1/;
   $line =~ s/\{/\[/g; $line =~ s/\}/\]/g;
   my $T = new Array<Set>(eval($line));
   my $edges = tropical_curve_from_subdivision($cayley_points, $T, $ptsindices);
   my $g = graph_from_edges($edges);
   my $cycle_length = 0;
   foreach my $comp (@{$g->BICONNECTED_COMPONENTS}) {if ($comp->size>2) {$cycle_length = $comp->size; last;}}
   ++$histogram{$cycle_length};
   my $ch = canonical_hash($g->ADJACENCY);
   print $OUT "ID :$id; edges: $edges; canonical_hash: $ch\n";
   if($id%10000 == 0){
      my $ttmp = Benchmark->new();
      my $diff = timestr(timediff($ttmp,$t0));
      print $LOG "$id($nclasses; $diff)\n";
   }
   if(defined $classes->{$ch}){
      my $test = graph::isomorphic($classes->{$ch}, $g->ADJACENCY) ? 1:0;
      # Small test to make sure there are no hash collisions.
      if($test == 0){
         die "Collision!\n";
      }
   } else {
      $classes->{$ch} = $g->ADJACENCY;
      $nclasses++;
   }
   ++$id;
}

$TDATA->close();
my $t1=Benchmark->new();

close $OUT or die "cannot close $output_file\n";

my $td1=timediff($t1,$t0);

print $LOG "\n$id processed lines\n", "elapsed time for computation: ", timestr($td1), "\n";
print $LOG "histogram:\n";
foreach my $os (sort {$a <=> $b} keys(%histogram)) {
    print $LOG $os, ":", $histogram{$os}, " ";
}
save($classes, $output_classes);
# if($nclasses != $classes->size()){ die "Something went wrong when collecting classes!\n"}
print $LOG "\nCollected $nclasses graphs\n";
close $LOG or die "cannot close $log_file\n";
